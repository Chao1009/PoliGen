<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# D7 — the grid structure-function backends become the shipped defaults: design

**Decision (author, 2026-10-27, `DECISIONS.md` D7):** the unpolarized default
moves from `ToyF2` to CT18NLO through LHAPDF (`--unpol-sf ct18nlo`), the
polarized default from `ToyG1` to NNPDFpol1.1 (`--pol-sf nnpdfpol`); the toys
stay as opt-in options; the leading-order evaluation of the two NLO sets is
filed as the development item. **This file is the design only.** Nothing
outside it was edited; no constant, tolerance, registry label or reference
number was moved. Every measurement below was made on 2026-10-27 on the tree at
`ca9b23a` (C++ build of 2026-10-26, `build/python/lipolgen/*.py` byte-equal to
`python/lipolgen/*.py`), LHAPDF 6.5.5 with `CT18NLO`, `NNPDFpol11_100` and
`EPPS21nlo_CT18Anlo_Li6`, PYTHIA 8.317, on 4 cores shared with three other
workers (load average 0.0–0.3 for the timing runs, up to 5.6 for the harness
runs); scratch scripts and outputs are in
`/tmp/claude-1000/-home-cpeng-Projects-polli-PoliGen/7790c59b-ea07-4677-a36f-0183e2410ef4/scratchpad/d7/`
(not permanent). **READ** marks what a file says (quoted with `file:line`),
**MEASURED** what was run, **INFERRED** what follows from them, and
**NOT VERIFIED** what could not be checked here. Line numbers are those of
`ca9b23a`: every `file:line` here cites the tree read for this design, and the
implementation will move many of them. (The per-section `Anchors as of`
markers this record first carried were taken out on 2026-10-27: the records
gate's as-of census, `test_every_as_of_anchor_passes_the_bar_at_its_commit`,
admits no new historical anchor, so every citation here is now checked
against the live tree, and one whose line the implementation moves fails that
gate until it is re-pointed (or annotated, which needs that census raised).)

Contents: §1 what the new defaults compute, and the letter's chain · §2 where
the default lives (the design) · §3 measured costs and behaviour changes ·
§4 inventory · §5 pin decisions · §6 benchmark rows · §7 published numbers ·
§8 CI and the missing-set failure mode · §9 implementation steps and gates ·
§10 verification plan · §11 risks · §12 questions only the author can answer ·
§13 the D7 decision text and the development item · §14 not verified.

---

## 1. What `--unpol-sf ct18nlo` and `--pol-sf nnpdfpol` compute

### 1.1 The unpolarized backend (READ)

* **The object.** `--unpol-sf ct18nlo` is `LhapdfSF("CT18NLO", 0)`, built one
  layer up by the binding (`python/bindings.cpp:4381-4390`, `set_unpol_sf`)
  because the core links no LHAPDF (`include/lipolgen/sf.hpp:200-201`,
  `include/lipolgen/pipeline.hpp:645-652`). It is **nucleon-level CT18NLO**:
  `F2p = x Σ e_q² (q + q̄)` over **five massless flavours** d u s c b with the
  proton charges {1/9, 4/9, 1/9, 4/9, 1/9}, the neutron by the u↔d charge swap
  (`src/lhapdf/lhapdf_sf.cpp:41-46`, `:59-69`) — the leading-order parton
  formula applied to an NLO set, no coefficient functions, no heavy-quark
  threshold, `max(F2, 0)` and `F2 = 0` outside 0 < x < 1 (`src/lhapdf/lhapdf_sf.cpp:61-68`).
  `f2n_over_f2p(x)` is evaluated at a fixed Q² = 10 (`src/lhapdf/lhapdf_sf.cpp:114-124`); it is what
  the T1 and T2 struck-nucleon species draws read.
* **The nucleus.** `NuclearF2::f2a = Z F2p + N F2n`, times an `emc_ratio(x)`
  hook only if one is given (`src/core/sf.cpp:127-131`). **No pipeline path
  sets that hook**: `default_inclusive_kernel` fills `opt.f2_source` and
  `opt.g1_model` and nothing else (`src/core/pipeline.cpp:1183-1188`), and the
  struck-cluster kernel forwards only `f2_source` / `g1_model`
  (`src/core/pipeline.cpp:283-288`). So a ⁶Li run on `ct18nlo` is
  3 F2p + 3 F2n of the **free** nucleon; ⁷Li is 3 F2p + 4 F2n; the tagged
  channels' struck clusters (`DEUTERON`, `TRITON`, `NEUTRON_TARGET`) the same
  free-nucleon sums; the coherent channel owns no structure function and
  rides the inclusive cell cross sections (`include/lipolgen/pipeline.hpp:606-617`).
* **EPPS21 is not applied, by any option.** `Epps21Ratio`
  (`include/lipolgen/lhapdf_sf.hpp:49-78`, `src/lhapdf/lhapdf_sf.cpp:147-190`)
  is a standalone class reachable from C++ and Python; no `PipelineConfig`
  field, CLI flag or kernel option reads it (grep of `include/ src/ python/`,
  2026-10-27). The one EPPS21 number inside a run is the polarized-EMC
  transfer scale `EMC_VALENCE_DEPLETION_EPPS21` (`include/lipolgen/sf.hpp:336`),
  which no kernel uses (`PolSF::g1_nucleus`'s `medium_ratio` is passed by
  nobody, `include/lipolgen/pipeline.hpp:700-701`).
* **F1, F_L and R.** `F1 = F2/(2x(1 + R))`, `F_L = F2 R/(1 + R)`
  (`include/lipolgen/sf.hpp:91-100`, `f1_from_f2` / `fl_from_f2`; `NuclearF2::f1a`, `src/core/sf.cpp:133-136`)
  with R = `r_sigma_lt` = 0.18/(1 + Q²/50) unless an `RFunc` is set
  (`include/lipolgen/sf.hpp:43-50`). The pipeline sets one only on
  `--b1-model li6-convolution --r-source r1998`
  (`include/lipolgen/pipeline.hpp:538-547`); **on the shipped `miller` path R
  stays the toy R under D7.** D7 moves F2; it does not move R.
* **Outside the grid.** CT18NLO's grid is x ∈ [1e−9, 1], Q ∈ [1.295, 1e5] GeV
  (`CT18NLO.info`: XMin, QMin, QMax), and LHAPDF's global `Extrapolator:
  continuation` (`deps/install/share/LHAPDF/lhapdf.conf`) **continues the
  evolution below Q² = 1.677 instead of freezing** (`include/lipolgen/sf.hpp:210-219`).
  The run measures it on every run: `meta["unpol_sf_grid_q2_min"]` and
  `meta["unpol_sf_below_grid_frac"]` (`python/bindings.cpp:686-693`, from
  `unpol_sf_grid_report`, `python/bindings.cpp:373-421`, against the channel's own rate), and
  prints it at the banner (`python/lipolgen/cli.py:1104-1124`). MEASURED
  (cell-centre rate weights, the classification `unpol_sf_grid_report` itself
  uses, `python/bindings.cpp:410-418`, `cell_rate_weights_pb`): 36.18 % of
  the inclusive ⁶Li config-1 rate below 1.677 (40.46 % at config 0, 34.17 % at
  config 2; ⁷Li config 1 36.23 %; tagged d–p 36.43 %; the coherent channel
  44.75 %, `include/lipolgen/pipeline.hpp:627-645`). No accepted cell centre
  is below x = 1.15e−4 or above Q² = 1897 GeV² at any shipped config — inside
  the grid in x and Q² everywhere else.

### 1.2 The polarized backend (READ)

* **The object.** `--pol-sf nnpdfpol` is `LhapdfG1("NNPDFpol11_100", 0)`
  (`python/bindings.cpp:4409-4419`): `g1 = ½ Σ e_q² (Δq + Δq̄)/x` over **three
  flavours** d u s (`src/lhapdf/lhapdf_sf.cpp:47-48`, `:71-82` `g1_scalar`) — leading order,
  no ΔC_q ⊗ Δq, no ΔC_g ⊗ Δg; Δc = Δb = 0 is NNPDFpol1.1's own convention
  (`src/lhapdf/lhapdf_sf.cpp:9-12`). The set is mem 0 = the average over 100 replicas
  (`NNPDFpol11_100.info`: SetDesc).
* **The nucleus.** `g1A = Z P_p g1p + N P_n g1n` with the per-nucleon
  effective polarizations of the ion (`src/core/sf.cpp:168-174`, `PolSF::g1_nucleus`), divided by A
  in the kernel (`src/core/xsec.cpp:88-90`, `InclusiveKernel::g1_nucleus`); g2 is Wandzura–Wilczek on the same
  g1 (`include/lipolgen/pipeline.hpp:698-699`). The A_∥ denominator F1 is the
  kernel's own (§1.1) — **CT18NLO with `r_sigma_lt`** under the D7 pair.
* **When the backends are mixed.** With `--pol-sf toy` the kernel builds
  `ToyG1` on its **own** `UnpolSF` (`src/core/xsec.cpp:80-82`), i.e. g1 =
  A1_toy(x) × F1(CT18NLO) under `--unpol-sf ct18nlo` — "`--pol-sf toy` does NOT
  mean g1 unchanged" (`include/lipolgen/pipeline.hpp:619-625`). With
  `--pol-sf nnpdfpol`, g1 is independent of the F1 backend; A_∥ = D g1/F1 is
  not.
* **Reach.** g1 is read only where the channel evaluates it (not
  `CoherentLi6`) **and** the fill carries λ_e·P_e ≠ 0
  (`src/core/pipeline.cpp:1337-1351`, `pol_sf_is_read`); elsewhere
  `meta["pol_sf"]` carries the label and `meta["pol_sf_reach"]` the reason
  (`python/bindings.cpp:675-684`). **At the CLI's default plan
  (`tensor-thirds`) `--pol-sf` is read on no channel**, so D7's polarized half
  changes nothing at the all-default run (MEASURED: `meta["pol_sf"]` = "not
  read under this run's fill").
* **Outside the grid.** NNPDFpol11_100 covers x ∈ [1e−5, 1], Q ∈ [1, 316.228]
  GeV (`NNPDFpol11_100.info`). MEASURED: no accepted cell centre of the
  shipped windows (⁶Li configs 0–2, ⁷Li, d) lies below Q² = 1 or x = 1e−5, so
  the polarized grid is not extrapolated at any cell centre (whether a
  per-event evaluation reaches below its cell centre's Q² was not checked).
  **There is no metadata for it** (`LhapdfG1` has no `q2_min`
  accessor; `unpol_sf_grid_report` handles `LhapdfSF` only,
  `python/bindings.cpp:403-421`): a user-lowered `Scenario::q2_min` would
  extrapolate g1 silently. §12 O-3.
* **Positivity.** The LO pair violates |A1| ≤ 1 at the window's top cell:
  A1 = g1/F1 = 2.02 (deuteron, x = 0.955, Q² = 1007) and 4.41 (free neutron,
  Q² = 1119) — MEASURED (§3.3); with the tensor term on the inclusive
  deuteron, that refuses six helicity combinations; the sampler refuses
  rather than clamp (`python/lipolgen/cli.py:1406-1432`).

### 1.3 Against the letter's `--pdf grid` chain (READ, sibling repository)

`--pdf grid` on the letter's four truth-level drivers means "F₂ᴬ from the
nuclear set EPPS21nlo_CT18Anlo_Li6 (not Z·F₂ᵖ + N·F₂ⁿ on a free proton), g₁
from NNPDFpol11_100, F₂ᵖ from CT18NLO for the g₁ denominator, and the SLAC/E143
R1998 fit in place of the toy R" (`PolarizedLithiumSim/plans/07_plb_letter_gluonometry.md:263-267`);
the code: `evgen/scripts/money_cos2phi.py:197-219` (R1998 on `grid`),
`fastsim/polli_fastsim/inputs.py:18-27`, `fastsim/polli_fastsim/structure.py:436-528`
(`NuclearF2FromGrid`), `fastsim/polli_fastsim/polarized.py:184-230`
(`PartonG1`). The EPPS21 grid in LHAPDF is the **isoscalar-averaged nucleon**
of ⁶Li (`docs/benchmarking/03_data_nuclear.md:107-111`: R_u = 0.578, R_d = 2.829
at x = 0.65 are the n/p average), so the sibling's proton-charge sum × A is
the whole-nucleus F₂ᴬ correctly.

| | LiPolGen under D7 (`ct18nlo` + `nnpdfpol`) | letter `--pdf grid` |
|---|---|---|
| F₂ᴬ of ⁶Li | 3 F₂ᵖ + 3 F₂ⁿ of free-nucleon CT18NLO, LO formula, 5 massless flavours | 6 × F₂ of EPPS21nlo_CT18Anlo_Li6 (nuclear: shadowing, EMC; on CT18ANLO), same LO formula and flavours (`structure.py:478`, :504-508) |
| below the grid's Q₀ | LHAPDF continuation below Q² = 1.677 (36.18 % of the 6Li config-1 rate), recorded in `meta` | F₂ᴬ **frozen** at Q₀² = 1.69 (`structure.py:502`); "8.7 % of the accepted cells but 36.3 % of the accepted one-year rate" (`plans/07:292-294`) |
| R (F₁, F_L, D(y)) | `r_sigma_lt` (toy R) | R1998 on every consumer |
| g₁ | `LhapdfG1`, NNPDFpol11_100, LO, d u s | `PartonG1`, the same set and formula |
| A_∥'s F₁ | CT18NLO with `r_sigma_lt` | EPPS21 with R1998 |
| ⁷Li | 3 F₂ᵖ + 4 F₂ⁿ, CT18NLO | no nuclear set: `NUCLEAR_F2_SET_FOR_ION` knows (3, 6) only (`structure.py:418`), `NuclearF2FromGrid` raises (:482-489) |
| tagged channels | the struck cluster on CT18NLO / NNPDFpol | the tagged drivers take no `--pdf` (toy) |
| PDF reader | LHAPDF 6.5.5 C++ | the pure-Python `parton` package, same grid files |
| measured σ ratio | 6Li config 1 ×0.7985 | N_DIS ×0.766 (`plans/07:352-357`, "within 4 %") |

**Not equal.** The D7 defaults are the letter's grid chain **without** the
nuclear set, **without** R1998 and **with** a continuation instead of a
freeze below Q₀. They share only g₁ and the LO F₂ formula.

**A finding about the premise, recorded, not acted on.** The decision's
wording "the grid backends the letter's production uses" (`DECISIONS.md` D7)
is wider than the letter's own text: its **published footing is the toy** —
"FOOTING. Single-fill estimator, toy structure functions with a placeholder R
… The two-fill estimator … and the nuclear-grid backend (plans/07 WP1
addendum) are bounded in Sec. 6 and are author calls, not adopted here"
(`paper/main.tex:28-31`; also `plans/07:361-372`, :1036-1042); the grid
column is run and quoted in §6 as a bound (`paper/main.tex:582-594`). The
`polligen` prototype's kernel defaults to `ToyF2`/`ToyG1`
(`evgen/polligen/xsec.py:454-469`). §13's entry text words this accurately;
the decision itself is the author's and is not reopened.

---

## 2. Where the default lives — the design

The core library cannot construct an LHAPDF object (`include/lipolgen/sf.hpp:200-201`).
Two places could carry the new default:

* **A (recommended) — the C++ `PipelineConfig` names the grids; the LHAPDF
  tier builds them.** `PipelineConfig::unpol_sf = UnpolSfSource::Ct18Nlo`,
  `pol_sf = PolSfSource::NnpdfPol` (`include/lipolgen/pipeline.hpp:986`,
  `include/lipolgen/pipeline.hpp:1003`); the objects stay empty in the core and are attached by the
  LHAPDF tier (new, §9 step 1) or by the binding. One definition of the
  default, read everywhere it is needed: `Pipeline::knob_provenance` already
  takes `at_default` from `const PipelineConfig d;` "the shipped defaults, read
  not retyped" (`src/core/pipeline.cpp:2117`, `:2678` `c.unpol_sf == d.unpol_sf`,
  `src/core/pipeline.cpp:2690` `c.pol_sf == d.pol_sf`), and the CLI and
  `make_config` read the default name from it instead of retyping `"toy"`.
  A C++ program that does not link the tier gets a refusal naming the tier,
  the set and `UnpolSfSource::Toy` — never the toy silently.
* **B (rejected) — the run surfaces default to the grids, the C++ struct stays
  `Toy`.** No C++ churn, but two definitions of one default (the struct and
  the CLI/Python), a `knob_provenance` `at_default` that is false on the
  shipped run unless special-cased, and a C++ `PipelineConfig{}` that runs the
  toy where Python runs CT18NLO — the entry-point dependence D7 is meant to
  end.

**Refinements of A** (each needed by a measured consequence in §3):

1. **Class-level null defaults stay the toy, and say so.**
   `InclusiveKernel::Options::f2_source/g1_model` (`include/lipolgen/xsec.hpp:198-199`,
   `src/core/xsec.cpp:59-87`), `StruckClusterOptions::f2_source/g1_model`
   (`include/lipolgen/pipeline.hpp:230-249`), `PythiaBridgeOptions::f2_source`
   (`include/lipolgen/pythia_bridge.hpp:280-284`), `BreakupOptions::f2`
   (`include/lipolgen/breakup.hpp:236-238`, null = the library's `ToyF2`), `ToyG1`'s base, the one-argument
   `default_inclusive_kernel(ion)` (`include/lipolgen/pipeline.hpp:1012-1017`)
   and `DeuteronConvolutionB1`'s default Options are building blocks the
   `Pipeline` always fills; their null branch is the toy because the core
   cannot be anything else. They keep it (the polligen port gate and
   `validation/reference/*.json` are pinned on exactly them, §5); their
   comments change from "the DEFAULT … every published number" to "the
   class-level toy; the run's default is set by `PipelineConfig`". The comment
   at `include/lipolgen/pipeline.hpp:1012-1013` ("The kernel the INCLUSIVE
   channel uses when `PipelineConfig::kernel` is null") becomes false under D7
   and must be reworded.
2. **One builder per named backend, in the LHAPDF tier.** New in
   `include/lipolgen/lhapdf_sf.hpp` / `src/lhapdf/lhapdf_sf.cpp`:
   `make_lhapdf_unpol_sf(UnpolSfSource)`, `make_lhapdf_pol_sf(PolSfSource)` and
   `attach_shipped_sf_backends(PipelineConfig&)` (fills every empty slot the
   tier can build: `unpol_sf_obj`, `pol_sf_obj`, and `b1_unpol_sf` on
   `Ct18Nlo`). The bindings' `set_unpol_sf` / `set_pol_sf` / `set_b1_unpol`
   (`python/bindings.cpp:4158-4196`, `:4364-4431`) call them, so "ct18nlo means
   `LhapdfSF("CT18NLO", 0)`" is written once (today it is written at
   `python/bindings.cpp:4177`, `:4383`, `:4411` (its `LhapdfG1("NNPDFpol11_100", 0)`) and in three `validate()`
   messages, `src/core/pipeline.cpp:1043-1078`).
3. **A pre-flight before any LHAPDF `Config` access** (§3.4, §8): every
   constructor in the tier checks `LHAPDF::findFile("lhapdf.conf")` and
   `LHAPDF::findFile("<set>/<set>.info")` **before** `lhapdf_quiet()`
   (`src/lhapdf/lhapdf_sf.cpp:86-89`, which calls `LHAPDF::setVerbosity` and
   so `Config::get()`), and throws one message naming the set, the searched
   paths, and the escape (`--unpol-sf toy` / `--pol-sf toy`, Python
   `unpol_sf="toy"`). `findFile` reads only the environment and the install
   prefix (`deps/src/LHAPDF-6.5.5/src/Paths.cc:19-53`).
4. **Python attaches; C++ asks.** The binding's `PipelineConfig()`
   constructor attaches the shipped objects **best effort** (never throws;
   a failed pre-flight leaves the slot empty and `validate()` reports the
   reason), so the documented `lg.PipelineConfig()` example
   (`python/lipolgen/__init__.py:9-15`) keeps working; `make_config(unpol_sf=None)`
   attaches **strictly** (throws the pre-flight message at once; the CLI
   turns it into a one-line `SystemExit`, `python/lipolgen/cli.py:1312-1313`).
   `cli.DEFAULTS["unpol_sf"]` / `["pol_sf"]` (`python/lipolgen/cli.py:823`)
   become `_l.unpol_sf_name(_l.PipelineConfig().unpol_sf)` /
   `_l.pol_sf_name(...)`, not a retyped string. Assigning
   `config.unpol_sf_obj = None` today resets to `Toy`
   (`python/bindings.cpp:4620-4628`); recommended under D7: it re-selects the
   shipped default (rebuilt by the tier), and the toy is chosen only by name.
5. **A process-wide cache of the two shipped objects** (optional): one
   `LhapdfSF("CT18NLO", 0)` and one `LhapdfG1("NNPDFpol11_100", 0)` per
   process, built under a mutex. It takes `make_config` from 30–36 ms back to
   ≈ 0.6 ms after the first call (§3.1). Sharing one LHAPDF object across
   generation threads is already what one `Pipeline` does; MEASURED bit-identical
   at 1 vs 4 threads (§3.6). Adopt only if the pytest wall time grows by more
   than 10 % without it.
6. **`validate()` under the new default** (`src/core/pipeline.cpp:1010-1145`):
   (a) a named backend with an empty slot (`src/core/pipeline.cpp:1028-1054`, `:1065-1081`): the
   message names D7, the set, `attach_shipped_sf_backends` and the toy escape;
   (b) **a caller-supplied `kernel`** (`src/core/pipeline.cpp:1104-1114`) is refused today beside
   any non-`Toy` selector — under D7 that refuses `make_config()` + `cfg.kernel
   = k` outright (MEASURED, §3.7). New rule: a kernel is accepted beside
   selectors **at their shipped default** (compared with `PipelineConfig{}`),
   the named-but-empty check is skipped then (the selectors are not read;
   `meta` already writes "caller-supplied kernel", `python/bindings.cpp:658-660`),
   and a non-default selector beside a kernel is refused as now;
   (c) **a directly set `struck.f2_source` / `struck.g1_model`**
   (`src/core/pipeline.cpp:1089-1102`) is refused today only when the selector is `Toy`; under any
   other selector the object silently wins while `meta` names the selector —
   MEASURED today (§3.5). New rule: refused unless the selector is `Custom`;
   (d) the `b1_unpol` collision (`src/core/pipeline.cpp:1132-1144`) stays; its message says that
   under the shipped default `--b1-model li6-convolution` needs `--b1-unpol`
   named (§12 O-1).
7. **The T2 bridge check moves into C++.** The species-draw consistency
   refusal lives only in the Python binding (`python/bindings.cpp:5107-5140`);
   a C++ caller binding a default-option `PythiaBridge` (ToyF2,
   `include/lipolgen/pythia_bridge.hpp:280-284`) to a D7 config draws T2
   species from the toy with no message (`examples/generate_full.cpp:152-173`
   does exactly that). Recommended: a C++ `set_pythia_hadronizer` in the
   PYTHIA tier that the binding calls.

---

## 3. Measured costs and behaviour changes (MEASURED, 2026-10-27)

### 3.1 Speed and set loading

Warm page cache (the set files had been read on 2026-10-26; a cold first read
was not measured — dropping the cache needs root). Fresh process per row.

| what | toy | `ct18nlo` + `nnpdfpol` |
|---|---|---|
| `import lipolgen` | 0.067–0.077 s | same |
| construct `LhapdfSF("CT18NLO", 0)`, first / a second one | — | 22.6–28.1 ms / 10.4–10.8 ms |
| construct `LhapdfG1("NNPDFpol11_100", 0)` | — | 12.6–13.8 ms |
| `make_config(events=10)` | 0.60–0.65 ms | 30.4–36.5 ms |
| `Pipeline` build, inclusive ⁶Li config 1 (3051 cells) | 0.055–0.059 s | 0.206–0.238 s |
| `Pipeline` build, tagged ⁶Li α / tagged ⁷Li α / coherent | 0.061–0.067 / 0.061 / 0.054 s | 0.221–0.227 / 0.215 / 0.208 s |
| T0 per event, inclusive ⁶Li, `tensor-thirds` (5000 ev, 2 runs) | 0.89 / 1.13 µs | 2.17 / 2.18 µs |
| T0 per event, inclusive ⁶Li, `helicity-flip` (5000 ev) | 0.94 µs | 2.40 µs |
| T0+T1 per event, tagged ⁶Li α, `tensor-thirds` (3000 ev, 2 runs) | 1.95 / 1.99 µs | 3.61 / 5.46 µs |
| T0+T1 per event, tagged ⁷Li α, `helicity-flip` (3000 ev) | 2.41 µs | 4.05 µs |
| per event, coherent ⁶Li (3000 ev) | 2.39 µs | 2.16 µs |
| `lipolgen-run --events 5000 --quiet`, end to end (3 runs) | 0.170–0.189 s | 0.352–0.398 s |

Reading: the switch costs ≈ 0.03 s of set loading and ≈ 0.15–0.18 s of
sampler set-up per run, and about ×2 per T0/T1 event (the species draws call
`f2n_over_f2p`, 20 `xfxQ2` per call); coherent is unchanged. At the CLI's
100 000-event default that is ≈ +0.3 s per run; the T2 (PYTHIA) tier (~40 k ev/s,
`docs/USAGE.md:2792`) dominates any hadronized run. `python/tests/test_throughput.py`'s
collapse floors (5e4 / 3e4 ev/s) pass at the new rates (simulated, §3.7).
Script: `.../scratchpad/d7/load_time.py`, `gen_speed.py`, `gen_speed.jsonl`.

### 3.2 One line of the price table, confirmed

`phase_A_baseline.md` §4.2 was measured at N = 200 000; the quantities are
N-independent (`../run_2026-10-26/phase_A_baseline.md:260-267`), so N = 2000 at seed 20261026
must reproduce them to the printed digits. It does
(`.../scratchpad/d7/confirm_price.py`):

| quantity | toy | `ct18nlo` + `nnpdfpol` | ratio | record |
|---|---|---|---|---|
| inclusive ⁶Li σ, `tensor-thirds` [pb] | 591846.1689 | 472571.9179 | ×0.79847 | 591846.1689 → 472571.9179 (×0.79847) |
| run-level A_zz | −5.193231e−4 | −6.503970e−4 | ×1.25239 | ×1.25239 |
| inclusive ⁷Li A_∥, `helicity-flip`, `--pol-sf nnpdfpol` alone | 5.291881e−4 | 2.812909e−4 | ×0.53155 | ×0.53155 |
| inclusive ⁷Li A_∥, **the D7 pair** | 5.291881e−4 | 3.534091e−4 | **×0.66783** | ×0.66783 |
| inclusive ⁷Li σ, `helicity-flip` | 590952.4264 | 470360.0619 | ×0.79594 | ×0.79594 |

The D7 price on ⁷Li A_∥ is the **pair's ×0.668**, not the ×0.53 of
`nnpdfpol` alone that `SUMMARY.md` and `DECISIONS_PRICED.md` §2 headline.
The run-level A_zz ratio is the inverse of the σ ratio to 6 × 10⁻¹²
(1/0.7984708573984506 = 1.2523938610085 against 1.2523938610006; MSTW the
same): on the shipped `miller` path the tensor shift does not depend on F₂,
so ×1.25239 on A_zz is the rate's normalisation and nothing else (INFERRED
from the two measured ratios).
Other rates measured on the way (default plans, config 1): coherent ⁶Li σ
12476.025 → 8806.096 pb (×0.70584); tagged d–p (`tensor-thirds`) 585044.23 →
456479.48 pb (×0.78025); tagged ⁷Li α (`helicity-flip`) 589760.77 → 467410.92 pb
(×0.79254); inclusive ⁶Li (`helicity-flip`) 591783.25 → 472509.00 pb
(×0.79845). The cos 2φ amplitude per cell does not move (F₂ cancels at fixed
Δ/F₁: 8.9e−16 over 3051 cells, `phase_A_baseline.md` §4.2), the super-bins by
≤ 0.23 %.

### 3.3 What the default pair refuses at the shipped window

The full matrix (`.../scratchpad/d7/matrix.py`): ⁶Li, ⁷Li and d × beam
configs 0, 1, 2 × every channel of the isotope × the four plans (⁷Li:
`helicity-flip` only, the spin-1 plans being refused at J = 3/2 by design) —
66 combinations, each built as a `Pipeline` at `x_max` = 1.0. **On the toy all
66 build; on the D7 pair 6 are refused by the sampler's positivity check:**

| isotope, config, channel, plan | refused cell | 1 + w_avg | A1 = g1/F1 |
|---|---|---|---|
| d, 0, inclusive, `helicity-flip` | x 0.955, Q² 283.7 | −0.004833 | 2.02 |
| d, 1, inclusive, `helicity-flip` | x 0.955, Q² 1007 | −0.01952 | 2.024 |
| d, 2, inclusive, `helicity-flip` | x 0.955, Q² 1897 | −0.02576 | 2.033 |
| d, 0, `tagged-d-p`, `helicity-flip` | x 0.955, Q² 229.7 | −0.005726 | 4.36 |
| d, 1, `tagged-d-p`, `helicity-flip` | x 0.955, Q² 1119 | −0.02414 | 4.406 |
| ⁶Li, 0, `tagged-6Li-alpha`, `helicity-flip` | x 0.955, Q² 480.9 | −0.04915 | 2.021 |

Only the d–p config-1 case was documented
(`python/lipolgen/cli.py:705-716`, `docs/USAGE.md:956-1004`); the inclusive-deuteron cases are new: inclusive d
`helicity-flip` is refused under `--unpol-sf ct18nlo` **alone** too
(1 + w_avg = −0.007475 at Q² = 1707), not under `--pol-sf nnpdfpol` alone.
`--x-max 0.95` cures both deuteron channels at config 1 (exit 0, 500 events);
the other four were not re-run at 0.95. The refused cell is the window's top
x = 0.955 cell in all six. On `tagged-d-p` and `tagged-6Li-alpha` (struck
kernels with no rank-2 slot by default) it is the vector term, where the LO
pair's A1 exceeds the positivity bound |A1| ≤ 1 (4.41, 2.02; §13, the
development item); on the inclusive deuteron the Miller b₁'s tensor term over
CT18NLO's small F₁ at x = 0.955 carries most of it — refused on `ct18nlo`
alone at A1 = 0.845 (INFERRED from the reported 1 + w_avg and A1; the split
was not computed). Further refusals the default now meets:

* **`--b1-model li6-convolution` without `--b1-unpol`**: refused by
  `validate()` (the `b1_unpol` collision, `src/core/pipeline.cpp:1132-1144`) —
  `--b1-unpol toy` means "the kernel's own UnpolSF", which is no longer ToyF2.
* **The mandatory band's ×2 row** (`--b1-band-scale 2`) on `cdks` and on
  `li6-convolution` (with `--b1-unpol ct18nlo` and with the gate's `mstw`) is
  refused at the documented `--x-max 0.95` **and** at 0.9 — cell x = 0.871,
  Q² = 186 (`cdks`, 1 + w_avg = −0.008234) / 315.3 (`li6-convolution`,
  −0.005282 on ct18nlo, −0.0109 on mstw) — and builds at **0.85**. Bands 0 and
  1 build at 0.95; every case builds on the toy at 0.95 (24 builds,
  `.../tasks/bjg0lxwat.output`). The recipe "RUN 0 / 1 / 2 and quote the
  envelope" at `--x-max 0.95` (`python/lipolgen/cli.py:534-546`,
  `docs/USAGE.md:619-641`) does not survive D7 as written.
* **`--inclusive-b1` on `tagged-6Li-alpha`**: refused at x_max 1.0
  (−0.1302 at x 0.955, Q² 167.3 — as documented, `python/lipolgen/cli.py:477-495`).

### 3.4 A missing set

Two environments, `LHAPDF_DATA_PATH='<dir>::'` (the trailing `::` blocks the
install-prefix fallback, `deps/src/LHAPDF-6.5.5/src/Paths.cc:27-31`):

* **`lhapdf.conf` present, sets absent:** `lipolgen-run --unpol-sf ct18nlo`
  prints `Info file not found for PDF set 'CT18NLO'` and exits 1;
  `--pol-sf nnpdfpol --plan helicity-flip` likewise for `NNPDFpol11_100`;
  `make_config(unpol_sf='ct18nlo')` raises `RuntimeError` with that text and
  the interpreter continues and exits 0. The set is named; the escape
  (`--unpol-sf toy`) and the searched path are not.
* **No `lhapdf.conf` (an empty data directory — the wheel's situation on any
  machine but the builder's, `docs/PACKAGING.md:334`):** the same commands
  **abort**: `terminate called after throwing an instance of
  'LHAPDF::ReadError' … Couldn't find required lhapdf.conf system config file`,
  SIGABRT, core dump, exit 134 — even when Python caught the first exception.
  Cause (READ): `Config::get()` constructs its function-static `_cfg` and then
  throws when `lhapdf.conf` is missing (`deps/src/LHAPDF-6.5.5/src/Config.cc:13-22`);
  at exit `~Config()` calls `verbosity()` → `Config::get()` → throws inside a
  destructor (:25-31, `deps/src/LHAPDF-6.5.5/include/LHAPDF/Config.h:56-58`) → `std::terminate`.
  LiPolGen reaches `Config::get()` through `lhapdf_quiet()` →
  `LHAPDF::setVerbosity(0)` in every tier constructor
  (`src/lhapdf/lhapdf_sf.cpp:39`, `:86-89`, `:98`, `:133`, `:157`). **Today
  only an explicit `--unpol-sf ct18nlo` / `--pol-sf nnpdfpol` can reach it;
  under D7 the all-default run does.** The toy default, which never touches
  LHAPDF, exits 0 in the same environment (MEASURED). The pre-flight of §2
  item 3 avoids `Config` entirely on this path.

### 3.5 A latent mislabel the default makes reachable

`make_config(channel="tagged-6Li-alpha", unpol_sf="ct18nlo")` with
`cfg.struck.f2_source = ToyF2()` set by hand is **accepted**:
`meta["unpol_sf"]` reads `"ct18nlo"` while σ = 591846.1688784088 pb, the toy
run's value to every digit (the struck kernel took the hand-set object,
`src/core/pipeline.cpp:1545-1546` fills only an empty slot, and `validate()`
guards the slot only under `Toy`, `src/core/pipeline.cpp:1089-1102`). Under D7 this needs no
`--unpol-sf` at all. §2 item 6(c) closes it.

### 3.6 Threads

On the D7 pair, `generate(0, False, 4)` is bit-identical to
`generate(0, False, 1)` in all 47 columns for `tagged-6Li-alpha`
(`tensor-thirds`), inclusive ⁶Li and tagged ⁷Li α (`helicity-flip`), 6000
events, seed 99 — one LHAPDF object shared by four generation threads. (A
test, not a proof of thread safety: `tests/test_b1_nuclear.cpp:150-152` says
an `LhapdfSF` is not to be shared across threads; LHAPDF 6.5.5's
`GridPDF` builds its interpolator lazily, `deps/src/LHAPDF-6.5.5/include/LHAPDF/GridPDF.h:210-214`,
and the `Pipeline` evaluates the kernel single-threaded at construction
before any generation thread starts.)

### 3.7 The Python suite under simulated D7 defaults

A scratch copy of the staged package (`.../scratchpad/d7/d7pkg/lipolgen/`, the
compiled `_lipolgen` symlinked) with two edits — `make_config` fills
`unpol_sf="ct18nlo"`, `pol_sf="nnpdfpol"` when given `None`, and
`cli.DEFAULTS` names them — run with `LIPOLGEN_TESTS_USE_INSTALLED=1`,
`PYTHONDONTWRITEBYTECODE=1`, `-p no:cacheprovider` (nothing written into the
repository) over the 22 test files that build runs. It simulates the **run
surfaces** of design A; it does not change the C++ struct (so a bare
`PipelineConfig()` stays toy in it) nor any harness. Baseline: the full suite
passed on this tree in phase F (`../run_2026-10-26/STATUS.md:348`, 1456 passed / 155 skipped / 0 failed).

| file | under simulated D7 | cause |
|---|---|---|
| `test_module.py` | 89 passed | — |
| `test_pipeline.py` | 18 passed | — |
| `test_cli.py` | 1 failed / 32 | `test_a_lower_case_isotope_still_runs[D-d]`: `--isotope D --plan helicity-flip`, refused at x 0.955 (§3.3) |
| `test_knob_provenance.py` | 17 failed / 490 / 155 skipped | 15 matrix cells vary `unpol_sf: toy -> ct18nlo` / `pol_sf: toy -> nnpdfpol` from a base that is now that value; 2 T2 cells build a `PythiaBridge` without `f2_source` (refused by `set_pythia_hadronizer`); 4 `li6-convolution` cells newly skip ("the BASE … is refused": `b1_unpol` toy under ct18nlo) |
| `test_sf_backend.py` | 27 failed / 27 | the file's gates are written on "the default is the toy" (G1, G2, G4, G7, G8, G9, G11, G12, G13) |
| `test_b1_model.py` | 15 failed, 4 errors / 32 | `li6-convolution` with `b1_unpol` toy under ct18nlo (14); a kernel beside a non-toy selector (1); the band fixture, whose `cdks` ×2 row at x_max 0.95 is refused (4 errors, §3.3) |
| `test_li7_rank2.py` | 3 failed / 24 | a caller kernel beside the default selectors (1); "the selectors reach the 7Li kernels" compares the default against `ct18nlo`, now equal (2) |
| `test_pzz_mode.py` | 3 failed / 22 | `PRICE` pins σ of run 2026-09-06 §A3 (591783.252… → 472509.001) |
| `test_rc_low_x_anchor.py` | 1 failed / 4 | the three-row band half-width table of run 2026-09-06 §A2 (max rel. diff 0.205) |
| `test_price_defaults.py` | 1 failed | the default row's σ, A_zz, A2 pins (the recipe prices against the default) |
| `test_pom_set_band.py` | 1 failed, 6 errors / 2 | the module's `PythiaBridge` without `f2_source` |
| `test_hfs.py` | 7 errors | the fixture's `PythiaBridge` without `f2_source` |
| `test_bench_chain_rc.py` | 2 failed / 2 | DJANGOH row history numbers (LiPolGen's own Born, §3.8) |
| `test_li6_unpol_band.py`, `test_meta_provenance.py`, `test_hepmc.py`, `test_cluster_config.py`, `test_release_metadata.py`, `test_throughput.py`, `test_rc.py`, `test_sampler.py`, `test_bench_t4_t6.py` | all pass (12, 8, 6, 27, 8, 5, 32, 9, 18 + 1 opt-in skipped) | — |

**71 failures and 17 errors in 9 files**; outputs in `.../scratchpad/d7/results/`.
Not covered by the simulation: the C++ suite (§4.4, §5.1 item 4), the bench-opt-in T4 row,
and the harness-test pins that follow §6's verdict switches.

### 3.8 Harness rows on the D7 Born / kernel

Run unchanged against the simulated package (`harness_*.{base,d7}.txt`):

* **`t4_djangoh_rad_noRad` (row 7), PASS → PASS.** The verdict R =
  σ_tail(LiPolGen)/σ_tail(DJANGOH) moves only through the tail table's
  interpolation of a ratio over a different Born: `polrad-full` 1.0199 /
  1.0261 / 1.0219 / 1.0098 → **1.0213 / 1.0262 / 1.0219 / 1.0090** (allowances
  0.0307 / 0.0310 / 0.0309 / 0.0308), `t-peak` 0.6118 / 0.8968 / 1.0112 /
  1.0067 → 0.6125 / 0.8969 / 1.0112 / 1.0068. The recorded sub-rows move:
  LiPolGen's own Born over DJANGOH's 1.469 / 1.116 / 1.155 / 1.065 → 1.124 /
  1.067 / 1.033 / 0.902, the ePIC-history ratios (t-peak 1.01617 → 1.02051 in
  the first bin; `polrad-full` 1.0271 → 1.03437, i.e. 47 % of DJANGOH's
  inelastic excess where `python/tests/test_bench_chain_rc.py:126-129` asserts
  "never more" than 40 %).
* **`t4_pythia_ep_closure`, PASS → PASS**: worst |bridge/stock − 1| 0.076
  (W2 n_ch_eta[−2, 0)) → 0.065 (W1 n_ch_eta[−∞, −2)); ⟨n_ch⟩ W1 1.007 → 1.009,
  W2 0.997 → 0.994 (the (x, Q²)-matched ratio; N = 20 000 per window).
* **`t6_hepmc3_convention`, PASS → PASS**: 23/23 clauses on the same 8
  `lipolgen-run` files (the d–p file is `tensor-thirds`, not refused).
* **`t5_delta_moment_bound` (row 13), PASS → PASS**, with the harness's
  `shipped_kernel()` replaced by the kernel a D7 `Pipeline` builds
  (`.../scratchpad/d7/delta_bound_d7.py`): `toy_1e-2[shipped]` saturation
  0.0644 / 0.0803 / 0.0941 → **0.0760 / 0.0907 / 0.1023** (V6 0.152 / 0.181 /
  0.205); helicity positivity w_min 0.9494 → 0.9135, A_T ∈ [−0.101, 0.051] →
  [−0.173, 0.087]. The letter's `moment_A` solved on the D7 F₁ gives A =
  −0.24038, **not** the printed −0.2924, so `a_reproduced` turns false: the
  letter models must keep the toy F₁, the letter's footing (§6).
* **`t5_weighted_unweighted`**: its REPORT line is unchanged in the
  simulation because it builds its own kernel with the one-argument
  `default_inclusive_kernel` (the toy); its header calls that "what `Pipeline`
  builds" (`validation/benchmarks/t5_weighted_unweighted.py:9-14`), false under
  D7. Not measured on the D7 kernel (the binding cannot pass backends to
  `default_inclusive_kernel`, `python/bindings.cpp:4720-4731`).

---

## 4. Inventory — every site that encodes the toys as the default

Found by `grep` over `include/ src/ python/ tests/ validation/ docs/ README.md
.github/ examples/ CMakeLists.txt` for `Toy`, `toy`, `UnpolSfSource`,
`PolSfSource`, `unpol_sf`, `pol_sf`, `unpol-sf`, `pol-sf` and the `default`
sentences around them (2026-10-27, `ca9b23a`), then read site by site. The
dated run records under `docs/open_items/run_*/` are history and are **not**
edited (the citation gate's `--records` pass keeps their links honest).
"toy" also names two things D7 does **not** touch — `toy_delta_gluon` (the
shipped Δ at 10⁻²) and `toy_b1` / `B1Mode::kToy` — and those sites are left
alone. Actions: **FLIP** the default · **REWORD** a "the toy is the default"
statement (no number) · **RETAKE** numbers with a dated D7 comment ·
**TOY** pass the toy explicitly with a dated D7 comment · **NEW** · **FIX** a
behaviour the switch makes reachable · **KEEP** (class-level toy or dated
history; reason given).

### 4.1 C++ (core and tiers)

| file:line | what | action |
|---|---|---|
| `include/lipolgen/pipeline.hpp:986` | `UnpolSfSource unpol_sf = UnpolSfSource::Toy;` | **FLIP** → `Ct18Nlo` |
| `include/lipolgen/pipeline.hpp:1003` | `PolSfSource pol_sf = PolSfSource::Toy;` | **FLIP** → `NnpdfPol` |
| `include/lipolgen/pipeline.hpp:977-1006` | the two fields' comments ("`Toy` is the DEFAULT and is today bit for bit") | REWORD (D7, the attach helper, the toy escape) |
| `include/lipolgen/pipeline.hpp:560-653` | `UnpolSfSource` block: "Toy the DEFAULT and bit for bit what every published number was made with" (579-583), the measured price, the grid clause | REWORD; the prices stay as toy-vs-grid prices |
| `include/lipolgen/pipeline.hpp:655-702` | `PolSfSource` block: "Toy the DEFAULT and bit for bit" (661-662) | REWORD |
| `include/lipolgen/pipeline.hpp:420-497` | `B1UnpolSource` block: "AT ITS SHIPPED DEFAULT `toy` that F1 is `ToyF2`'s" (467-468), the gate's "DEFAULT `ToyF2`" (403) | REWORD (`--b1-unpol` keeps `toy` as its own default, §12 O-1) |
| `include/lipolgen/pipeline.hpp:946-964` | `b1_unpol` field comment | REWORD |
| `include/lipolgen/pipeline.hpp:230-249`, `:255-258` | `StruckClusterOptions::f2_source/g1_model` "BOTH NULL BY DEFAULT … `ToyF2`"; `struck_cluster_kernel` "on the library's default toy F2 and toy g1" | KEEP the null = class-level toy; REWORD "default" → "class-level"; say the `Pipeline` fills them at the shipped D7 defaults |
| `include/lipolgen/pipeline.hpp:1012-1017` | one-argument `default_inclusive_kernel(ion)`: "The kernel the INCLUSIVE channel uses when `PipelineConfig::kernel` is null" | KEEP the toy (port gate, `b1_default_li6.json`); REWORD — false under D7 |
| `include/lipolgen/pipeline.hpp:1030-1047` | the 8-argument form: "NULL — the default … `ToyF2`" | REWORD (null = class-level toy; the `Pipeline` passes the D7 objects) |
| `include/lipolgen/pipeline.hpp:724-740` | `PipelineConfig::kernel` comment (selectors vs a caller kernel) | REWORD for the new rule (§2 item 6b) |
| `src/core/pipeline.cpp:1019-1081` | `validate()` named-but-empty / toy-with-object clauses and their messages | **FIX**: D7 message; skip beside a caller kernel at default selectors |
| `src/core/pipeline.cpp:1089-1102` | the struck-slot rule, `Toy` only | **FIX** → refuse unless `Custom` (§3.5) |
| `src/core/pipeline.cpp:1104-1114` | kernel beside a non-`Toy` selector refused | **FIX** → compare with `PipelineConfig{}` |
| `src/core/pipeline.cpp:1132-1144` | the `b1_unpol` collision | KEEP the rule; REWORD the message (D7) |
| `src/core/pipeline.cpp:1163-1188` | `default_inclusive_kernel` comment: "NULL -- the default … keeps the shared `ToyF2`" | REWORD |
| `src/core/pipeline.cpp:2117`, `:2665-2690` | `knob_provenance`: `at_default` from `PipelineConfig d` ("the shipped defaults, read not retyped"), compared in the `unpol_sf` / `pol_sf` rows | KEEP — becomes right by itself once the struct flips |
| `src/core/pipeline.cpp:1373-1424` | `pol_sf_reach_report`: "at pol_sf = toy the kernel's own ToyG1 is built on the kernel's own UnpolSF" | KEEP (true of the explicit toy) |
| `src/core/pipeline.cpp:2055-2072`, `:2600` | dated toy σ in comments (12476.025113 → 12446.109030 pb) | KEEP (measured statements with their configuration); add "(toy)" where they say "default" |
| `include/lipolgen/sf.hpp:12-13`, `:18-20` | "a *toy* implementation always available"; R "the toy default" | KEEP (true; R is not D7) |
| `include/lipolgen/sf.hpp:210-219` | `LhapdfSF::q2_min` doc: "a `--unpol-sf ct18nlo` run puts a third of its rate there" | REWORD → "the shipped default run" |
| `include/lipolgen/sf.hpp:226-238` | `LhapdfG1` | **NEW** `q2_min()` accessor (§12 O-3) |
| `include/lipolgen/xsec.hpp:178-181`, `:198-199`; `src/core/xsec.cpp:59-85` | `Options::f2_source` "default: ToyF2", `g1_model` "default: ToyG1 on f2_source" | KEEP (class level); REWORD "default" → "null" |
| `include/lipolgen/pythia_bridge.hpp:280-284` | "Null = the library's default `ToyF2`" | KEEP; REWORD; **NEW** C++ `set_pythia_hadronizer` check (§2 item 7) |
| `include/lipolgen/breakup.hpp:236-238` | `BreakupOptions::f2` null = `ToyF2` | KEEP (the `Pipeline` fills it from the kernel) |
| `include/lipolgen/b1_nuclear.hpp:27-41` | the gate's "SHIPPED DEFAULT unpolarised backend, `ToyF2`, gives 0.440" | REWORD (it is `--b1-unpol`'s default; under D7 `li6-convolution` must name it) |
| `include/lipolgen/lhapdf_sf.hpp:34-39`, `src/lhapdf/lhapdf_sf.cpp:86-100`, `:132-135`, `:154-162` | `lhapdf_quiet()` before any check; `mkPDF` | **FIX** pre-flight (§2 item 3); **NEW** builders + `attach_shipped_sf_backends` |
| `src/lhapdf/lhapdf_sf.cpp:1-23` | the file header ("Ports … of the LHAPDF-grid classes") | REWORD: add "the shipped defaults since D7; LO evaluation — the development item" |
| `CMakeLists.txt:82-110` | LHAPDF not found → `message(WARNING …)` and the tier silently off | REWORD the warning: "the shipped defaults (D7) need this tier; a build without it refuses every default run" (§8) |

### 4.2 Bindings, Python package, CLI

| file:line | what | action |
|---|---|---|
| `python/bindings.cpp:4198-4230`, `:4233-4270` | enum docstrings "Toy the DEFAULT …" | REWORD |
| `python/bindings.cpp:4364-4431`, `:4158-4196` | `set_unpol_sf` / `set_pol_sf` / `set_b1_unpol` build `LhapdfSF`/`LhapdfG1` inline | route through the tier builders; messages name the escape |
| `python/bindings.cpp:4523-4549` | `StruckClusterOptions` property docstrings "None (the default) is the kernel class's own ToyF2" | REWORD (class level) |
| `python/bindings.cpp:4604-4655` | `PipelineConfig.unpol_sf` / `pol_sf` docstrings ("Toy is the default"); `unpol_sf_obj` / `pol_sf_obj` setters (None → Toy, `:4620-4650`) | REWORD; **FIX** None → the shipped default (§2 item 4) |
| `python/bindings.cpp:4552-4553` (`py::class_<PipelineConfig>`, `py::init<>()`) | `py::init<>()` | **NEW** best-effort attach (§2 item 4) |
| `python/bindings.cpp:4720-4731` | `default_inclusive_kernel` binding passes `nullptr, nullptr` for the two backends | **NEW** `unpol_sf=None, pol_sf=None` arguments (harnesses need the shipped kernel) |
| `python/bindings.cpp:403-421`, `:686-693` | `unpol_sf_grid_report`, `meta["unpol_sf_below_grid_frac"]` | KEEP (0 on the toy, 0.3618 on the shipped 6Li run); **NEW** the polarized twin if O-3 is taken |
| `python/bindings.cpp:640-660` | meta comments ("toy and ct18nlo differ by a FACTOR 0.80") | KEEP |
| `python/bindings.cpp:5107-5140` | `set_pythia_hadronizer` species-draw check | move into C++ (§2 item 7) |
| `python/lipolgen/__init__.py:248-287` | `UNPOL_SF`: "\"toy\" is the DEFAULT and is bit for bit" | REWORD |
| `python/lipolgen/__init__.py:289-323` | `POL_SF`: "\"toy\" is the DEFAULT and is bit for bit" | REWORD |
| `python/lipolgen/__init__.py:174-203` | `B1_UNPOL`: "at its shipped default \"toy\" it is `ToyF2`" | REWORD |
| `python/lipolgen/__init__.py:497-526`, `:668-684` | `make_config` doc ("Both default to the toy backends and are then bit for bit") and the selector code (`_lipolgen.set_unpol_sf` / `set_pol_sf`; None leaves the C++ default) | REWORD; **FLIP** None → attach the C++ default strictly |
| `python/lipolgen/__init__.py:828-844` | `run()`: "None (the `unpol_sf=\"toy\"` default) leaves the bridge's own ToyF2" | REWORD (the bridge now gets CT18NLO by default) |
| `python/lipolgen/__init__.py:8-16` | module docstring builds `lg.PipelineConfig()` and runs it | KEEP (works with §2 item 4) |
| `python/lipolgen/cli.py:823` | `DEFAULTS`: `unpol_sf="toy", pol_sf="toy"` | **FLIP**, read from `_l.PipelineConfig()` |
| `python/lipolgen/cli.py:614-650` | `--unpol-sf` help: "'toy' (default) is the library's ToyF2" (622) | REWORD (generated from the C++ default name) |
| `python/lipolgen/cli.py:651-688` | `--pol-sf` help: "'toy' (default) is ToyG1" (667) | REWORD |
| `python/lipolgen/cli.py:547-576` | `--b1-unpol` help: "at its shipped default 'toy' the SPIN-BLIND cell cross section is bit for bit" (564-566) | REWORD; add the D7 refusal |
| `python/lipolgen/cli.py:477-495`, `:534-546`, `:695-725` | `--inclusive-b1`, `--b1-band-scale`, `--x-max` helps: "OFF THE TOY UNPOLARISED BACKEND", the d–p case, "RUN 0 / 1 / 2" at 0.95 | REWORD with §3.3's measured cases (deuteron, ⁶Li config 0, band ×2 at 0.85) |
| `python/lipolgen/cli.py:877-920` | `B1_UNPOL_GATE_ROW`, `b1_gate_lines`: "The SHIPPED DEFAULT toy F2 is 0.440" (902-903) | REWORD ("the toy F2, `--b1-unpol toy`") |
| `python/lipolgen/cli.py:1004-1047` | `require_unpol_sf_tier` / `require_pol_sf_tier`: tier only | **NEW** the set pre-flight's message (§8) |
| `python/lipolgen/cli.py:1049-1125` | `sf_banner_lines`: "BOTH ARE THE TOY BACKENDS" (1088-1097) | KEEP for an explicit toy run; **NEW** a default-run clause naming the LO caveat (the development item) |
| `python/lipolgen/cli.py:1379-1388` | the bridge's `f2_source = cfg.unpol_sf_obj`: "bit for bit at `--unpol-sf toy`" | KEEP the code; REWORD |
| `python/lipolgen/cli.py:1406-1432` | `Pipeline` refusal comment (the d–p case) | REWORD (six default-path cases) |
| `python/README.md:212-213`, `:262-266` | "a toy backend is always there"; the bridge note | REWORD (D7 default) |

### 4.3 Examples

| file:line | what | action |
|---|---|---|
| `examples/generate_tagged.cpp:139`, `examples/generate_coherent.cpp:125`, `examples/generate_full.cpp:152` | `PipelineConfig cfg;` then a `Pipeline` | **NEW** `#include "lipolgen/lhapdf_sf.hpp"` + `attach_shipped_sf_backends(cfg)` (CMake skips them when the tier is off, `CMakeLists.txt:186-200`) |
| `examples/generate_full.cpp:169-173` | `PythiaBridge(bc)` built with default options (no `f2_source`) | **FIX** `bopt.f2_source = cfg.unpol_sf_obj` |
| `examples/generate_inclusive.cpp:100-110` | its own `InclusiveKernel::Options` (class-level toy) | attach the D7 objects (it is the C++ API page); the `docs/USAGE.md` §8 throughput table it measures is re-measured or labelled "toy" |

### 4.4 C++ tests (`tests/`)

Every C++ test that builds a `Pipeline` default-constructs a
`PipelineConfig`; none names a backend (`grep unpol_sf|pol_sf tests/*.cpp`:
only `Custom` objects in `test_pipeline.cpp:1349-1547` and `b1_unpol` in
`test_b1_nuclear.cpp`). Under design A every one of them meets "a named
backend with an empty slot". Root construction sites (copies of a root
inherit it):

| file | sites | what they test | action |
|---|---|---|---|
| `tests/test_pipeline.cpp` | 33 `PipelineConfig`, 25 `Pipeline` | machinery (conservation, closure, threading, share invariance, HepMC3, hooks, optics); the port gate vs the Python prototype (`:398-437` RP tag fractions 0.0247 / 0.3061 / 0.2545 at 10 %, `:438-472` the 7Li optics-blind fractions); the D2 cases (`:1377-1547`, incl. `:1467-1469` `CHECK_NOTHROW(cfg.validate()); // the shipped default`) | **TOY** via one helper for machinery and port gate; **RETAKE** the D2 cases for D7 (the default refuses unattached; `Custom` cases unchanged); **NEW** LHAPDF-guarded D7 cases (§9 step 2) |
| `tests/test_b1_nuclear.cpp` | 17 `PipelineConfig`, 37 `validate()` | T16 `--b1-unpol` reach and refusals (`:2053-2213`), T18 `--r-source` (`:2341-…`); G3b on the library `ToyF2` (0.440) | **TOY** (they test the b1 selectors and pin the toy row the gate is read against); **NEW** one D7 case: `li6-convolution` under the shipped default with `b1_unpol` toy is refused with the D7 message |
| `tests/test_rc_pipeline.cpp` | 7, via `base_cfg` (`:41`) | RC band/tail through the pipeline | **TOY** in `base_cfg` (one site) |
| `tests/test_t2.cpp` | 10 | the T2 bridge on pipeline events | **TOY** (the bridge's own null `f2_source` is the toy, so the pair stays consistent) |
| `tests/test_breakup.cpp` | 6, via `tagged_cfg` (`:40`) | T1 breakup | **TOY** in `tagged_cfg` |
| `tests/test_fsi.cpp` | 5 (`PipelineConfig cfg;` at `:385`, `:511`, `:605`) | FSI weight | **TOY** |
| `tests/test_triton_sf.cpp` | 2 (`PipelineConfig cfg;` at `:177`) | Ciofi–Simula through the pipeline | **TOY** |
| `tests/test_reference.cpp`, `test_sf.cpp`, `test_xsec.cpp`, `test_tagged.cpp`, … | kernels/samplers built directly (class-level null backends) | the polligen port gate (`validation/reference/*.json`) | KEEP (class-level toy unchanged under A) |
| `tests/test_lhapdf.cpp` | the LHAPDF classes directly | `LhapdfSF`/`LhapdfG1`/`Epps21Ratio` against `parton` | KEEP; **NEW** pre-flight cases (missing set → the D7 message, no `Config` access) |

The helper: one header `tests/sf_test_support.hpp` with
`inline PipelineConfig toy_config()` (sets `unpol_sf = Toy`, `pol_sf = Toy`)
and a dated comment — "a machinery test; the shipped D7 default needs the
LHAPDF tier's objects, which the core suite does not build; D7's own wiring
is gated in `test_lhapdf.cpp` / the D7 cases of `test_pipeline.cpp`" — so the
replacement `PipelineConfig cfg;` → `PipelineConfig cfg = toy_config();` is
mechanical and greppable.

### 4.5 Python tests (`python/tests/`) — the per-test decisions are §5

| file:line | what | action |
|---|---|---|
| `test_sf_backend.py` (whole file; header `python/tests/test_sf_backend.py:1-56`) | the D2 gates written on "the default is the toy" | **RETAKE** as D2+D7 (§5.2) |
| `python/tests/test_knob_provenance.py:298-301`, `:312-315` | `V("unpol_sf", "toy -> ct18nlo")`, `V("pol_sf", "toy -> nnpdfpol")`, the kernel variant, the `b1_unpol` cells on `li6-convolution` | **RETAKE** the two SF variants (`ct18nlo -> toy`, `nnpdfpol -> toy`); FOLLOW the `CONV` cells (`CONV` names `b1_unpol`) |
| `test_knob_provenance.py` (the T2 cells) | `PythiaBridge` without `f2_source` | FOLLOW: `opts.f2_source = cfg.unpol_sf_obj` |
| `test_price_defaults.py:39-50` | `SIGMA_TOY_PB`, `AZZ_TOY`, `A2_TOY`, `RATIOS` | **RETAKE** (the recipe's default row) |
| `test_cli.py:163-178` | `--isotope D --plan helicity-flip` | FOLLOW: add `--x-max 0.95` |
| `test_li7_rank2.py:210-230`, `:329-…` | caller kernel; "selectors reach the 7Li kernels" | KEEP under §2 item 6b / **RETAKE** (toy vs default) |
| `test_pzz_mode.py:168-195` | `PRICE` (run 2026-09-06 §A3) | **TOY** |
| `test_rc_low_x_anchor.py` (`test_the_three_row_table_of_band_half_widths_on_azz`) | run 2026-09-06 §A2 table | **TOY** |
| `python/tests/test_b1_model.py:66-80` (`_cfg`), `:684-727`, `:808-809` | the b1 selectors on the toy; `DEFAULTS["b1_unpol"] == "toy"` | **TOY** in `_cfg`/`_unpol_cfg`; the `b1_unpol` default asserts stay (unchanged default) |
| `test_pom_set_band.py:52` (`_run`) | hand-built bridge; README's measured DPDF band | **TOY** + wire the bridge |
| `test_hfs.py` (fixture) | hand-built bridge; HFS schema | wire the bridge (follow the default) |
| `test_bench_chain_rc.py:107-175` | DJANGOH row history numbers and the "< 0.4" bound | **RETAKE** with the harness (§6, recommended) |
| `test_bench_t1_nucleon.py:272-442` | `configs["toy"]`, `headline == "toy"`, `status == configs["toy"]["status"] == "fail"` | **RETAKE** with the harnesses' verdict switch |
| `python/tests/test_bench_t5_delta_bound.py:341-369`, `:458`, `:497-508`, `:598-601` | `toy_1e-2[shipped]` saturation, positivity, REPORT text; `a_reproduced is True` | **RETAKE** the shipped row; the letter rows stay on the toy kernel (a_reproduced stays true) |
| `python/tests/test_bench_t5_closures.py:100-101`, `:236`, `:264` | `bjorken[toyg1]` as the shipped sub-row (`subrows`; `bjorken[nnpdfpol11]` beside it) | **RETAKE** the label |
| `test_bench_t1_r_world.py:323-562` | `f2_backends["toy"]` labels/status | **RETAKE** labels only (verdict is R's) |
| `test_bench_t4_t6.py` (bench-opt-in T4 pins) | worst 0.076 | **RETAKE** (0.065 measured, §3.8) |
| `test_bench_report.py` | the committed `REPORT.*` `--check` | follows the regeneration |
| `conftest.py:90-103` (`polligen` fixture) | the sibling prototype | KEEP (only `test_hfs.py`'s schema reader uses it; the numeric port gate is `validation/reference/`) |
| `conftest.py` | — | **NEW** session check: the LHAPDF tier is on and the shipped sets are found, else one `pytest.exit` naming them (§8) |

### 4.6 Harnesses and validation scripts

| file:line | what | action |
|---|---|---|
| `validation/benchmarks/t1_nmc_f2d.py:21-23`, `:72-73`, `:167`, `:344`, `:439`, `:487` | "toy ToyF2() -- the SHIPPED default … THE REPORT ROW'S VERDICT IS THIS CONFIGURATION'S"; `HEADLINE = "toy"`; labels | **HEADLINE → `ct18nlo`**; labels; dated note (§6) |
| `validation/benchmarks/t1_nmc_f2d_over_f2p.py:33-35`, `:81-82`, `:150`, `:353`, `:452`, `:505` | the same rule ("ToyF2 (shipped default)", `HEADLINE = "toy"`) | **HEADLINE → `ct18nlo`** (FAIL → PASS) |
| `validation/benchmarks/t1_hermes_f2d.py:19-21`, `:170`, `:353`, `:515`, `:606` | the same rule, on F₂ᵈ ("ToyF2 (shipped default)", `HEADLINE = "toy"`) | **HEADLINE → `ct18nlo`** |
| `validation/benchmarks/t1_g1d_world.py:77-78`, `:165`, `:235`, `:473`, `:594` | "ToyG1 on ToyF2 -- the SHIPPED default … THE REPORT ROW'S VERDICT IS THIS ONE'S" | **HEADLINE → `nnpdfpol`** |
| `validation/benchmarks/t5_sum_rules.py:8`, `:64-78`, `:184` | "Bjorken sum on the shipped `ToyG1`"; "the failing sub-row is the shipped toy's" | REWORD (the shipped sub-row becomes `bjorken[nnpdfpol11]`); `bc_ww` stays on `ToyG1()` (an implementation closure) |
| `validation/benchmarks/t5_delta_moment_bound.py:66-73`, `:628-630`, `:865` | "SHIPPED (the verdict): `default_inclusive_kernel(li6())`'s Delta … F1 the shipped ToyF2" | `shipped_kernel()` → the D7 kernel; **NEW** `letter_kernel()` = the toy, for `moment_A`/`moment_B` (they reproduce the letter's printed A = −0.2924 on its toy footing) |
| `validation/benchmarks/t1_r_world.py:52-53`, `:446`, `:605` | "the SHIPPED ToyF2's F2p" (a recorded H1 variant); `k_ship = default_inclusive_kernel(li6())` | REWORD labels (the toy variant stays recorded, a CT18NLO one is "shipped"); `k_ship` → the D7 kernel (the a₂ ratio is F₂-free) |
| `validation/benchmarks/t4_djangoh_rad_noRad.py:20-21`, `:33-35`, `:269-271`, `:590-592` (`lg.make_config`) | "LiPolGen's own Born is `InclusiveKernel::dsigma_unpol` with the shipped `UnpolSfSource.Toy`" — **outside** the hashed block (`:41-187`, `DESIGN_BEGIN`/`DESIGN_END` `:369-370`) | recommended: follow the shipped default — re-word `:20-21` and `:33-35`, add a dated D7 addendum after `:187` (`DESIGN BLOCK END`; the block is never edited), re-take the recorded sub-rows (§6); alternative: `make_config(..., unpol_sf="toy", pol_sf="toy")` at `:590` |
| `validation/benchmarks/t4_pythia_ep_closure.py:60-61` | "(… LiPolGen's kernel on ToyF2)" | REWORD (CT18NLO since D7); the row follows the default |
| `validation/benchmarks/t5_weighted_unweighted.py:9-14`, `:143-145` | "the shipped defaults … `default_inclusive_kernel(li6())` (Miller b1, ToyF2, ToyG1), which is what `Pipeline` builds" | the D7 kernel (binding extension, §4.2); RETAKE |
| `validation/benchmarks/t2_hermes_b1_table2.py:15-17`, `:34`, `:160` | `cdks_conv_default` = the class's default Options (`ToyF2`); F₂ for the inversion = `_l.MstwSF()`, `_l.ToyF2()` if no PYTHIA | KEEP (by design; not the run default) |
| `validation/benchmarks/t3_nmc_li6_over_d.py` | `Epps21Ratio` on CT18NLO | KEEP |
| `validation/price_defaults.py:521`, `:736`, `:769`, `:817`, `:864`, `:1068-1074` | the default row "unpol_sf = toy (shipped default)"; family defaults | **RETAKE**: the default row = the D7 pair (and D5/D6's defaults), `unpol_sf = toy` / `pol_sf = toy` become priced rows |
| `validation/dump_polligen_reference.py:12`, `validation/dump_b1_default_li6.py` | dumps on the toy backends / the one-argument kernel | KEEP (the port gate; the class-level toy) |
| `validation/vmc_tag_fractions.py:140-147` | `make_config` defaults | follows the default; the RP tag fraction carries no structure function (INFERRED; not re-run) |

### 4.7 Documents

| file:line | what | action |
|---|---|---|
| `README.md:39-41` | "the default run is the one every reference JSON is pinned against at rtol 1e-12" | REWORD (the reference JSONs pin the class-level toy kernels) |
| `README.md:56-58` | selector table: `--b1-unpol` "against 0.440 on the shipped `toy`"; `--unpol-sf` / `--pol-sf` "`toy`, bit for bit" | **FLIP** the default column; the measured prices stay |
| `README.md:160-167`, `:186-198` | "the toy costs ×0.80 on σ(⁶Li)…"; "on the default toy F₂ it is 0.440" | REWORD |
| `README.md:218-222`, `:243-256` | dependencies; the wheel "the pure-C++ generator works … `LhapdfSF` fails" | REWORD: the default needs the three sets; the toy run is what works without data |
| `docs/USAGE.md:11-67` | pip install, data paths | add: the shipped defaults need `CT18NLO`, `NNPDFpol11_100` on `LHAPDF_DATA_PATH`, else run `--unpol-sf toy --pol-sf toy` |
| `docs/USAGE.md:405-420`, `:536-618`, `:619-641`, `:642-742` | §2a gate and b1 tables: "The shipped default unpolarised backend — the library `ToyF2`" (409), "the default `toy`" (417, 548, 567), "Every row was made with the default `ToyF2`" (720), the top-x-cell recipe | REWORD; add §3.3's band ×2 at 0.85 |
| `docs/USAGE.md:743-781` | §2b: "Both default to the toy backends and are then bit for bit" (758-759) | REWORD |
| `docs/USAGE.md:782-881` | §2b "What the toy has been costing, measured" | KEEP as the toy-vs-grid price; retitle |
| `docs/USAGE.md:935-1005`, `:1056-1085` | the d–p `--x-max` subsection; "What is refused" table | REWORD for §3.3 and §2 item 6 |
| `docs/USAGE.md:1140-1276` | §2c ⁷Li: σ = 590952.42641509 pb for both categories (1162) | label "(toy)" or re-measure on the default (470360.0619 pb, §3.2) |
| `docs/USAGE.md:2768-2793` | §8 throughput table (2.83 / 1.02 / 2.32 M ev/s) | re-measure with the D7 examples or label "toy" |
| `docs/PHYSICS_CHANNELS.md:120`, `:136-143`, `:191`, `:196`, `:198-199` | "default `ToyG1`"; "The production defaults — `ToyF2`, `ToyG1`, `r_sigma_lt`"; "both defaulting to the toy and then bit for bit" | REWORD (`r_sigma_lt` stays the default R); citation gate (§9) |
| `docs/PHYSICS_CHANNELS.md:128` | 7Li σ 590952.42641509 pb measured | label "(toy)" |
| `docs/OPEN_ITEMS_SOLUTIONS.md:68`, `:92` (row 14, `sf-backend-injection`), `:2592` | state board and §14 "OPT-IN SHIPPED … both default to the toys" | dated D7 note: "shipped as the default 2026-10-27" |
| `docs/benchmarking/08_channel_priority.md:327`, `:333`, `:361-364` | ingredient A ("unpolarized nucleon inputs F₂, R") and H rows' development-item column; §5 item 3 "The toy F₂ / g₁ failures (A/H) … are the decisions this file asks the author to take at M2" | dated D7 note: taken; the A/H development item = §13's text; §7 revision log |
| `docs/benchmarking/BENCHMARK_PLAN.md:98-120` | dated status 2026-09-26 ("FAIL recorded for the shipped ToyF2 / ToyG1") | KEEP; add a dated D7 status paragraph |
| `docs/benchmarking/00_in_tree_checks.md:244` | "shipped `ToyF2` default, 0.440" | dated note |
| `docs/benchmarking/02_data_nucleon_deuteron.md`, `docs/references/01_nucleon-deuteron-data.md`, `docs/references/REFERENCES.md` | 11 / 12 / 7 mentions of the shipped toy | dated note where "shipped" is claimed |
| `docs/PACKAGING.md:334` | "`LhapdfSF('CT18NLO')` FAILS … then LHAPDF's own `terminate`" | REWORD: now the default path; the pre-flight replaces the abort |
| `validation/benchmarks/README.md:224-235` | status cells (F-2 "FAIL, recorded") | F-2 → PASS on the shipped `ct18nlo`; the others' verdict labels |
| `validation/README.md:83-97` | the reference dumps on `ToyF2`/`ToyG1` | KEEP (class-level toy) |
| `docs/CONVENTIONS.md`, `docs/HEPMC3_CONVENTION.md` | no default-dependent example (only `pol_sf` labelling, `docs/CONVENTIONS.md:37`, `:59` (`--pol-sf` under the unpolarised-beam plans); the xsec description, `HEPMC3_CONVENTION.md:181-189`) | KEEP |

### 4.8 CI, packaging, environment

| file:line | what | action |
|---|---|---|
| `.github/workflows/ci.yml:73`, `.github/scripts/build_deps.sh:62`, :137-153 | the three sets fetched from `/current/`, **unpinned** (`docs/PACKAGING.md:55-66`) | **NEW** sha256 pins for the three tarballs (§8) |
| `.github/workflows/ci.yml:423-438` | wheel smoke step `env -i … lg.run(channel='inclusive', isotope='6Li', events=20000, seed=1)` | **RETAKE** (§8): toy run with nothing set; default run refused cleanly under `LHAPDF_DATA_PATH=<empty>::`; default run with the sets |
| `.github/workflows/ci.yml:240-252` | C++ suite, pytest, citation and SPDX gates | KEEP; **NEW** a "sets absent" step (§8) |
| `.github/workflows/ci.yml:332-339` | PYTHIA-tier-off job asserts `HAVE_LHAPDF` | KEEP |
| `env.sh` | `LHAPDF_DATA_PATH=$LIPOLGEN_DEPS/share/LHAPDF` | KEEP |

---

## 5. Pin decisions

### 5.1 The rule

The run's ground rule (`STATUS.md:25-28`): a number that depends on the
decision is re-taken with a dated comment naming it, or the test passes the
old option explicitly when its purpose is to pin the old physics. Applied:

1. **RETAKE** — the pin's purpose is a property of the run surface's
   **default** (what the default run gives, records, prints or prices): the
   D2 selector gates, the knob-provenance cells, the price recipe's default
   row, the harness verdict rows whose declared rule follows "the shipped
   default". Comment: `# D7 (2026-10-27): re-taken on the shipped
   ct18nlo/nnpdfpol defaults; was <old> on ToyF2/ToyG1.`
2. **TOY** — the pin reproduces a number in a **dated record** or a
   **port-gate reference**, or pins the toy's own physics (G3b 0.440): pass
   `unpol_sf="toy", pol_sf="toy"` (C++ `UnpolSfSource::Toy`,
   `PolSfSource::Toy`). Comment: `# D7 (2026-10-27): the toy named explicitly
   -- this pins <record>, made on ToyF2/ToyG1.`
3. **FOLLOW** — a machinery, schema, identity or closure test with no
   backend-specific number runs on the default unchanged; if a D7 consequence
   refuses it (the top-cell positivity, the `b1_unpol` naming, the bridge
   wiring), make the smallest change that keeps its purpose (add
   `--x-max 0.95`, name `--b1-unpol`, wire `f2_source`), with the comment.
4. **C++ core suite** — rules 2 and 3 collapse into `toy_config()`, because
   the core suite builds no LHAPDF object; the D7 wiring gets its own
   LHAPDF-guarded cases (§9 step 2).

### 5.2 Per pin

| test | decision | the new value / the reason |
|---|---|---|
| `test_sf_backend.py` G1 "the default is bit for bit" | RETAKE | split: G1a the **toy** named explicitly is bit for bit the class-level null branch (the pre-D7 default, all 47 columns, three channels); G1b the default equals `ct18nlo`/`nnpdfpol` named explicitly, bit for bit |
| `test_sf_backend.py` G2–G4 (enum/object in step; named-but-empty; caller kernel; b1 collision; struck slot) | RETAKE | §2 items 4, 6a–6d: `unpol_sf_obj = None` → the shipped default; kernel accepted at default selectors, refused beside non-default; struck slot refused unless `Custom` (§3.5 as the regression case) |
| `test_sf_backend.py` G5 (CLI surface) | RETAKE | `DEFAULTS` = the C++ default names; help says `'ct18nlo' (default)` |
| `test_sf_backend.py` G6, G9, G12, G13 (the tagged kernel moves; pol_sf not on coherent; the d–p pair refusal; bridge on another backend) | RETAKE | compare **toy vs default** instead of default vs grid; G12's "each backend alone never needed it" re-measured on the default pair |
| `test_sf_backend.py` G7, G10, G11 (meta, grid fraction, kernel on tagged) | RETAKE | `meta["unpol_sf"] == "ct18nlo"` at default; below-grid 0.3618 / 0.4475; G11 unchanged in logic |
| `test_knob_provenance.py:298-301` | RETAKE | variants `ct18nlo -> toy`, `nnpdfpol -> toy`, `lhapdf=True` kept (the base needs the sets) |
| `test_knob_provenance.py:312-315` and the `CONV` cells | FOLLOW | `CONV` gains `b1_unpol="ct18nlo"`; the toy-path cell becomes `unpol_sf="toy", b1_unpol="toy"` |
| `test_knob_provenance.py` T2 cells | FOLLOW | wire `opts.f2_source = cfg.unpol_sf_obj` |
| `test_price_defaults.py:39-50` | RETAKE | default row σ 472571.9179123315 pb, A_zz −6.503970486e−4, A2(S1..S4) from the D7 run; `RATIOS` gains the toy as a priced row (σ ×1.2523939, A_zz ×0.7984709 against the new default) — all taken from the run itself, after I4–I6 if they precede (§9) |
| `test_cli.py:163-178` | FOLLOW | `--x-max 0.95` on the inclusive `helicity-flip` command (all spellings) |
| `test_li7_rank2.py:210-230` | KEEP | passes under §2 item 6b; else TOY |
| `test_li7_rank2.py:329-…` (both) | RETAKE | toy (explicit) vs the default: σ ×0.79594 / ×0.79254 |
| `test_pzz_mode.py:173-180` `PRICE` | TOY | run 2026-09-06 §A3's σ (591783.2520093301 …) |
| `test_rc_low_x_anchor.py` three-row table | TOY | run 2026-09-06 §A2's half-widths |
| `test_b1_model.py` (`_cfg`, `_unpol_cfg`, band fixture, gate rows) | TOY | the b1 selectors and G3b 0.440 are toy statements; **NEW** two D7 cases: `li6-convolution` alone is refused under the default; the ×2 band row needs `--x-max 0.85` on the default |
| `test_pom_set_band.py` | TOY + wire | README's DPDF band (−2.6 % / +9.9 %, kaons ×2.8) was measured on the toy |
| `test_hfs.py` fixture | FOLLOW | wire the bridge |
| `test_bench_chain_rc.py:107-175` | RETAKE (recommended) | §3.8's numbers; the "(v − 1)/(e − 1) < 0.4" bound restated as measured (0.474 in the first bin on the CT18NLO Born), with the reason |
| `test_bench_t1_nucleon.py:272-442` | RETAKE | `headline == "ct18nlo"` (F₂ᵈ, ratio) / `"nnpdfpol"` (g₁ᵈ); `dp["status"] == "pass"` |
| `test_bench_t5_delta_bound.py` | RETAKE (shipped row) + TOY (letter rows) | 0.07599861 / 0.09072656 / 0.10229499; w_min 0.9134945; `a_reproduced` stays true on the toy kernel |
| `test_bench_t5_closures.py` | RETAKE | the shipped Bjorken sub-row's label |
| `test_bench_t1_r_world.py` | RETAKE | labels only |
| `test_bench_t4_t6.py` (T4 opt-in pins) | RETAKE | worst 0.065 (§3.8; re-run at the implementation) |
| `tests/test_pipeline.cpp:398-472` (port gate vs the Python prototype) | TOY | the prototype is toy-based (`evgen/polligen/xsec.py:454-469`) |
| `tests/test_pipeline.cpp:1377-1547` (the four `TEST_CASE("D2: …")` blocks) | RETAKE | as `test_sf_backend.py` |
| every other C++ `PipelineConfig` root (§4.4) | TOY (`toy_config()`) | machinery; the core suite builds no LHAPDF object |
| `validation/reference/*.json` consumers (`tests/test_reference.cpp`, `python/tests/test_kernel_reference.py`, `test_module.py`) | KEEP | class-level null backends = the toy; optionally name `ToyF2()` in `_kernel_6li` for belt and braces |

**Which published numbers move** is §7.

---

## 6. Benchmark rows: declared rules and the expected flips

| row | declared rule about "the shipped backend" (READ) | today's verdict row | under D7 | status |
|---|---|---|---|---|
| `t1_nmc_f2d` | "toy ToyF2() -- the SHIPPED default … THE REPORT ROW'S VERDICT IS THIS CONFIGURATION'S" (`validation/benchmarks/t1_nmc_f2d.py:21-23`); "The row's STATUS is the shipped default's (toy)" (`:72-73`); `HEADLINE = "toy"` (`:167`) | toy 8307.9/155 FAIL | `ct18nlo` 3401.8/146 FAIL (toy recorded) | FAIL → FAIL |
| `t1_nmc_f2d_over_f2p` | the same ("the SHIPPED default": `validation/benchmarks/t1_nmc_f2d_over_f2p.py:33-35`, `:81-82`; `:150` `HEADLINE`) | toy 2161.4/211 FAIL | `ct18nlo` 215.5/189, p = 0.0904 | **FAIL → PASS** |
| `t1_hermes_f2d` (row 14) | the same, on F₂ᵈ ("the SHIPPED default": `validation/benchmarks/t1_hermes_f2d.py:19-21`, `:170`) | toy 5303.5/68 FAIL | `ct18nlo` 1788.2/45 FAIL | FAIL → FAIL |
| `t1_g1d_world` | "ToyG1 on ToyF2 -- the SHIPPED default … THE REPORT ROW'S VERDICT IS THIS ONE'S" (`validation/benchmarks/t1_g1d_world.py:77-78`, `:165`, `:235` `HEADLINE`) | toy 107.8/15, 69.4/15, 54.9/28 FAIL | `nnpdfpol` 45.0/15, 52.6/15, 51.0/28 FAIL | FAIL → FAIL |
| `t5_sum_rules` | "bjorken -- A MEASUREMENT OF THE SHIPPED MODEL … on `ToyG1` this sub-row measures what the shipped default g1 does" (`validation/benchmarks/t5_sum_rules.py:64-78`); "fail if any evaluated sub-row fails" (`:181-183`) | `bjorken[toyg1]` FAIL (divergent), `bjorken[nnpdfpol11]` FAIL | the shipped label moves to `bjorken[nnpdfpol11]` (+19.3 / +12.4 / +9.6 %, not converged) | FAIL → FAIL |
| `t5_delta_moment_bound` (row 13) | "SHIPPED (the verdict): `default_inclusive_kernel(li6())`'s Delta … F1 the shipped ToyF2" (`validation/benchmarks/t5_delta_moment_bound.py:66-73`) | 0.0644 / 0.0803 / 0.0941 PASS | 0.0760 / 0.0907 / 0.1023 (MEASURED, §3.8); positivity w_min 0.9135 | PASS → PASS |
| `t4_djangoh_rad_noRad` (row 7) | "LiPolGen's own Born is … the shipped `UnpolSfSource.Toy`" (`validation/benchmarks/t4_djangoh_rad_noRad.py:20-21`, outside the hashed block); "The verdict is on R … which needs no Born at all" (`:29-30`) | `polrad-full` 1.0199 / 1.0261 / 1.0219 / 1.0098 PASS | 1.0213 / 1.0262 / 1.0219 / 1.0090 (MEASURED) | PASS → PASS |
| `t1_r_world` (row 12) | "THE ROW'S VERDICT IS THIS BACKEND'S" — R's, `r_sigma_lt` (`validation/benchmarks/t1_r_world.py:11-18`); the ToyF2 F₂ᵖ is a RECORDED H1 variant (`:52-53`) | 12.54/10 PASS | unchanged (D7 does not move R; the a₂ ratio is F₂-free) | PASS → PASS |
| `t4_pythia_ep_closure` | "LiPolGen's kernel on ToyF2" (`validation/benchmarks/t4_pythia_ep_closure.py:60-61`); verdict at matched (x, Q²) | worst 0.076 PASS | worst 0.065 (MEASURED) | PASS → PASS |
| `t5_weighted_unweighted` | "the shipped defaults … `default_inclusive_kernel(li6())` (… ToyF2, ToyG1), which is what `Pipeline` builds" (`validation/benchmarks/t5_weighted_unweighted.py:9-14`) | pull +0.359 PASS | an identity; expected PASS on the D7 kernel | NOT VERIFIED |
| `t6_hepmc3_convention` | `lipolgen-run` with defaults | 23/23 PASS | 23/23 (MEASURED) | PASS → PASS |
| `t2_hermes_b1_table2` (row 3), `t3_nmc_li6_over_d` (row 4) | F₂ by design (`MstwSF`; `Epps21Ratio` on CT18NLO) | — | unchanged | — |

**The report: 13 pass / 10 fail → 14 / 9** (only `t1_nmc_f2d_over_f2p`
flips), if D4–D6 move nothing else. Every verdict switch is the harness
following its own declared rule, not a new rule; the toy stays a recorded
sub-row with its numbers. The headers get a dated D7 note outside any hashed
block, and the REPORT's plan-row texts follow `make_report.py`'s regeneration
(`validation/benchmarks/make_report.py:473` `plan_row_map`, `:648-654`). The status cells of
`validation/benchmarks/README.md:224-235` and the dated status paragraph of
`BENCHMARK_PLAN.md` follow.

---

## 7. Published numbers that move, and what regenerates them

| where | numbers | regenerated by |
|---|---|---|
| `docs/benchmarking/REPORT.md`, `REPORT.json`, `figures/` | §6 | `LIPOLGEN_BENCH=1 python3 validation/benchmarks/make_report.py --all --plots --compare <pre-I7 REPORT.json> --out-dir docs/benchmarking` |
| the price table (`../run_2026-10-26/phase_A_baseline.md` §4.2 is dated, stays) | a new default row and the toy as a priced row | `python3 validation/price_defaults.py --plot` → the I7 record |
| `README.md:56-58` selector table | the default column; the prices are toy-vs-grid prices and stay | by hand |
| `README.md:39-41` ("pinned against at rtol 1e-12"), `:160-167` ("the toy costs ×0.80"), `:186-198` ("0.440"), `:243-256` (the wheel's `LhapdfSF` failure) | sentences, no new number | by hand |
| `docs/USAGE.md:782-881` price tables ("What the toy has been costing, measured") | stay (prices, both ends named) | — |
| `docs/USAGE.md:1162`, `docs/PHYSICS_CHANNELS.md:128` | ⁷Li σ 590952.42641509 pb at "the default" | by hand: label "(toy)" or re-take 470360.0619 pb |
| `docs/USAGE.md:2785-2787` | throughput 2.83 / 1.02 / 2.32 M ev/s | the examples, re-run (machine-dependent) or label "toy" |
| `docs/USAGE.md:2585-2590` | RC rate cost 575 k → 495 k ev/s | dated; keep with "(toy)" |
| `docs/CONVENTIONS.md`, `docs/HEPMC3_CONVENTION.md` | none depend on the default (§4.7) | — |
| `validation/benchmarks/README.md:224-235`, `BENCHMARK_PLAN.md` | status cells, a dated paragraph | by hand, from the regenerated report |

No reference JSON moves (`git diff --stat validation/reference/` must stay
empty): every dump is a class-level toy object.

---

## 8. CI and the missing-set failure mode

**What CI has (READ).** `.github/workflows/ci.yml:73` `LHAPDF_SETS: "CT18NLO
NNPDFpol11_100 EPPS21nlo_CT18Anlo_Li6"`, part of the cache key (:79-80);
`.github/scripts/build_deps.sh:62`, :137-153 fetches them from
`lhapdfsets.web.cern.ch/current/` into the prefix, and the `test` job
restores that cache with `fail-on-cache-miss: true` (`ci.yml:181-190`) and
exports `LHAPDF_DATA_PATH` from `env.sh` (:208-223). So the sets are present
wherever the C++ suite and pytest run in CI — the suite already needed them
on 2026-09-26 (13 C++ and 39 pytest failures, every one "Info file not found
for PDF set", `../run_2026-09-26/SUMMARY.md:242-243`). The workflow has never
run on GitHub (`ci.yml:3-28`; registry row 14), so none of this is verified
there. The sets are **unpinned** ("`/current/` is a moving target by
construction", `docs/PACKAGING.md:61-66`).

**The rule.** Results must not depend on the environment: the default either
runs on CT18NLO / NNPDFpol11_100 or refuses by name. Never the toy in silence.

**The design.**

1. **The pre-flight** (§2 item 3), in every LHAPDF-tier constructor, before
   `lhapdf_quiet()`:
   `findFile("lhapdf.conf")` and `findFile("<set>/<set>.info")`; on a miss,
   `std::runtime_error` with, verbatim in substance: *"the shipped default
   `--unpol-sf ct18nlo` (D7, 2026-10-27) needs the LHAPDF set 'CT18NLO', and
   LHAPDF cannot find it (searched: <LHAPDF::paths()>; LHAPDF_DATA_PATH=<value
   or unset>; lhapdf.conf <found at … / NOT found>). Install it
   (`.github/scripts/build_deps.sh` fetches the three sets), point
   LHAPDF_DATA_PATH at it, or run the toy explicitly with `--unpol-sf toy`
   (Python `unpol_sf="toy"`) and say so — it is never substituted."* The same
   for `NNPDFpol11_100` / `--pol-sf toy`. This also removes the exit-time
   abort of §3.4 (no `Config` is touched on the refusal path). The CLI already
   turns a `RuntimeError` from `make_config` into a one-line `SystemExit`
   (`python/lipolgen/cli.py:1312-1313`); a set problem therefore exits 1 with
   that line, never a traceback or a core dump.
2. **The tier absent** (`-DLIPOLGEN_WITH_LHAPDF=OFF`, or `lhapdf-config` not
   found): `HAVE_LHAPDF` is false; `PipelineConfig()` leaves the named slots
   empty; every default run refuses with the tier message
   (`python/lipolgen/cli.py:1004-1047` / `src/core/pipeline.cpp:1028-1081`,
   reworded); the toy runs. `CMakeLists.txt:107-109` turns a missing
   `lhapdf-config` into a WARNING and builds on: under D7 that WARNING says the
   build cannot run the shipped defaults (a configure-time FATAL_ERROR unless
   `-DLIPOLGEN_WITH_LHAPDF=OFF` is passed explicitly is the stricter
   alternative; recommended, so no one ships a default-less build by accident).
3. **The C++ suite without the sets.** The core suite builds no LHAPDF object
   (§4.4, `toy_config()`), so it runs; the LHAPDF cases (`tests/test_lhapdf.cpp`,
   the new D7 cases) **fail** with the pre-flight's message — loud, as today's
   13, but one clear line each instead of LHAPDF's bare "Info file not found".
   A missing set on a tier-ON build is an environment fault, not a skip
   (skipping would let a suite pass that never ran the shipped default).
   Without the tier those cases are not compiled (`#ifdef LIPOLGEN_HAVE_LHAPDF`,
   `tests/test_lhapdf.cpp:34`).
4. **pytest without the sets.** Nearly every run-level test goes through the
   default, so a missing set would surface as hundreds of identical
   failures. `conftest.py` gains a `pytest_sessionstart` check: if
   `lg.HAVE_LHAPDF` and a shipped set fails the pre-flight (a new
   `_lipolgen.lhapdf_set_problem(name) -> str | None`, the same text), the
   session stops at once with `pytest.exit(<the message>, returncode=3)`.
   `LIPOLGEN_TESTS_ALLOW_MISSING_SETS=1` lets it proceed (each default-path
   test then fails with the pre-flight's message; toy-only tests pass) — for
   a developer who knows. With the tier OFF the same stop, naming the tier.
   CI never sets the variable.
5. **Pin the sets.** `build_deps.sh` verifies the tarballs it fetches; the
   artefacts this tree was measured against (sizes equal to
   `docs/PACKAGING.md:55-57`'s `content-length`):
   `CT18NLO.tar.gz` `c9127231e77e97cbec79cb5839203ab00f8db77237a061b61f9420f2b7b9c213`,
   `NNPDFpol11_100.tar.gz` `48b586b14821ef2858a52a763d3db99ff283a1405d7677298d4c9ec12f0dde5a`,
   `EPPS21nlo_CT18Anlo_Li6.tar.gz` `1afe5f6f8df61edb158485b419a788528f4c1eed1d9769e57eebc59c665cd868`
   (`sha256sum ../deps/src/*.tar.gz`, 2026-10-27). Optionally `meta` records
   the loaded member's identity (`unpol_sf_set`: name, member, `DataVersion`,
   sha256 of `CT18NLO_0000.dat` = `375db856…a72e5410`; `NNPDFpol11_100_0000.dat`
   = `c64ba5fe…c45f96`), so a run made on a moved upstream set says so.
6. **CI jobs.** (a) `test`: unchanged; the defaults now run on the restored
   sets. (b) **NEW step "shipped sets absent"** in `test` (cheap, no rebuild):
   with `LHAPDF_DATA_PATH=$RUNNER_TEMP/empty::` — once with a copied
   `lhapdf.conf`, once without — assert `lipolgen-run --events 10` exits 1 (not
   134) with "CT18NLO" and "--unpol-sf toy" in stderr, and
   `lipolgen-run --unpol-sf toy --pol-sf toy --events 10` exits 0; and pytest
   stops with return code 3 naming the sets. (c) **wheel**
   (`ci.yml:423-438`): the smoke `lg.run(...)` under `env -i` would hit the
   vendored `libLHAPDF`'s compiled-in prefix — which on the runner is the
   restored cache path, so it would pass for the wrong reason. Replace it by:
   `lg.run(..., unpol_sf='toy', pol_sf='toy')` with nothing set (the wheel runs
   without data); the default under `LHAPDF_DATA_PATH=$RUNNER_TEMP/empty::`
   raising the pre-flight `RuntimeError` and the interpreter exiting 0 after
   catching it; the default with `LHAPDF_DATA_PATH` pointing at the restored
   sets. (d) `pythia-tier-absent`: unchanged (it keeps the LHAPDF tier,
   `ci.yml:332-339`).
7. **Packaging.** The wheel vendors `libLHAPDF` but not the set store or
   `lhapdf.conf` (`docs/PACKAGING.md:334`, registry row 15): under D7 a
   pip-installed default run needs the user's LHAPDF data. `docs/USAGE.md`
   §"Installing the Python package with pip" and `docs/PACKAGING.md` say so,
   with the toy escape. Registry row 15 (publish the wheel or not) gains this
   consequence; deciding it stays the author's.

---

## 9. Ordered implementation steps (I7) and the gate after each

**Position.** I7 runs **after I4, I5 and I6**, sequentially, as its own
commit, so that every number it re-takes is attributable to D7 alone
(`include/lipolgen/pipeline.hpp:244-249` gives the reason in this tree's own
words: two independent movements shipped together are unattributable). D7
moves every rate and touches the same files as the others
(`test_knob_provenance.py`, the README selector table, `price_defaults.py`,
`t4_djangoh_rad_noRad.py`, which D5 relabels for the tail); going last means
its re-takes are final and the price table's default row is measured once on
the D4–D6 tree. Save the pre-I7 `docs/benchmarking/REPORT.json` before
step 1.

Gates (numbers from `docs/open_items/run_2026-10-26/STATUS.md:176-183`, `:346-353`):
**B** = `cmake --build build -j2` (zero warnings) + `./build/lipolgen_tests`
(≈ 7 min on 4 shared cores); **P** = `python3 -m pytest python/tests -q`
(≈ 8 min); **Pb** = `LIPOLGEN_BENCH=1 python3 -m pytest python/tests -q -k bench`
(≈ 5 min); **G** = `python3 validation/check_physics_channels_links.py`
(strict), `--records-only`, `--audit-ranges` — on a MOVED block or a drifted
named point, `--fix` re-points it; on an EDITED block, re-read the citing row
and only then `--record-ranges` (R5 refuses a block the sentence does not
evidence); never `--loose` as the gate (`validation/check_physics_channels_links.py:184-206`,
`validation/check_physics_channels_links.py:236-241`); **S** = `python3 validation/check_spdx_headers.py`;
**R** = `LIPOLGEN_BENCH=1 python3 validation/benchmarks/make_report.py --all --plots --compare <pre-I7 REPORT.json> --out-dir docs/benchmarking`,
then `make_report.py --check`.

1. **The LHAPDF tier.** Pre-flight; `make_lhapdf_unpol_sf` /
   `make_lhapdf_pol_sf` / `attach_shipped_sf_backends`; `LhapdfG1::q2_min`
   (if O-3); the optional cache. New `tests/test_lhapdf.cpp` cases: a missing
   set and a missing `lhapdf.conf` (via `LHAPDF::setPaths("<empty>::")`,
   restored after the case) give the message and touch no `Config`; the builders return the same objects the
   bindings built. A Python subprocess test: `LHAPDF_DATA_PATH=<empty>::`
   `make_config()` raises naming CT18NLO and `unpol_sf="toy"`, and the process
   exits 0. **Gate: B.**
2. **The core default and `validate()`.** FLIP `include/lipolgen/pipeline.hpp:986`, `:1003` (`UnpolSfSource::Toy`,
   `PolSfSource::Toy`);
   §2 item 6 rules; the comment blocks of §4.1. C++ tests: `toy_config()`
   everywhere §4.4 says TOY; the D2 cases re-taken; NEW LHAPDF-guarded cases —
   the default refuses unattached and runs attached, `meta` reads ct18nlo /
   nnpdfpol, a kernel at default selectors is accepted, a struck slot under a
   non-`Custom` selector is refused (§3.5's configuration), threaded
   generation on the grids is bit-identical (§3.6). **Gate: B** (every case
   passes; the tally changes by the new cases only; `git diff --stat
   validation/reference/` empty).
3. **Bindings, package, CLI.** Builders in `set_*`; `PipelineConfig()`
   best-effort attach; `unpol_sf_obj = None`; `default_inclusive_kernel`
   gains the two backends; `make_config` strict attach; `cli.DEFAULTS` from
   C++; helps, banners, refusal texts (§4.2); C++ `set_pythia_hadronizer`
   (§2 item 7); conftest's session check. Python tests per §5.2 (all but the
   bench ones). **Gate: B, P.**
4. **Harnesses and the price recipe.** §4.6: the four `HEADLINE`s, the
   Bjorken label, `t5_delta_moment_bound`'s two kernels, `t1_r_world`,
   `t4_djangoh_rad_noRad` (outside the hashed block; its `design_block_sha256()`
   must still match, `validation/benchmarks/t4_djangoh_rad_noRad.py:448-466`), `t4_pythia_ep_closure`,
   `t5_weighted_unweighted`, `price_defaults.py`; their bench tests (§5.2).
   **Gate: Pb, then R** — every `--compare` change must be in §6's table
   (14 / 9 if D4–D6 moved nothing else), and nothing else may move.
5. **Examples.** §4.3; run each once (`--events 2000`) and, if the §8
   throughput table is to be re-measured, then. **Gate: B.**
6. **Documents.** §4.7, one dated D7 note per site, the numbers from steps 2–5
   only. **Gate: G, S**, and `P` again for `test_doc_link_gate.py`.
7. **CI and packaging.** §8 items 5–7 (`build_deps.sh` pins; the two new CI
   steps; the wheel step). The workflow cannot be run here: run the
   "sets absent" commands locally with the same `LHAPDF_DATA_PATH` values.
   **Gate: S;** the local re-enactment.
8. **The record** `phase_I7_backends.md`: every re-taken number with old →
   new, the `--compare` output, the gate tallies, and the measured §3 table
   re-run on the implemented tree. **Gate: G** (`--records-only`).

---

## 10. Verification plan

* **Bit-identity of the toy.** Before and after I7, the same seeded run with
  `--unpol-sf toy --pol-sf toy` (inclusive ⁶Li `tensor-thirds`, tagged ⁶Li α,
  coherent, inclusive ⁷Li `helicity-flip`, 2000 events each): sha256 over all
  47 columns, `sigma_pb` and `sigma_per_category_pb` equal — the toy option is
  the old default, bit for bit.
* **The new default equals the named grids.** After I7, a default run and the
  same run with `ct18nlo`/`nnpdfpol` named: identical hashes, and identical to
  the pre-I7 tree's named-grid run (D4–D6 being in both).
* **The numbers of §3.2 reproduce** on the implemented tree at N = 2000 (σ
  472571.9179, A_zz −6.503970e−4, ⁷Li A_∥ 3.534091e−4, coherent 8806.096 pb),
  times whatever D4–D6 changed (D4 does not touch these; D6 does not reach the
  inclusive channel).
* **The six refusals of §3.3** happen at the default with the cell and the
  cure in the message, and `--x-max 0.95` runs all six; the band ×2 row runs
  at 0.85 and is refused at 0.95.
* **The missing-set behaviour** of §8 item 6(b), both environments, exit 1
  and never 134.
* **The reference JSONs**: `git diff --stat validation/reference/` empty.
* **The report**: `--compare` against the pre-I7 report shows exactly §6's
  changes; `--check` OK.
* **An adversarial pass** (the run's phase V) reads this plan's §4 against
  the diff: every listed site touched or explicitly kept, every re-taken pin
  carrying the dated comment, no tolerance or declared rule edited (diff the
  harness headers' rule paragraphs and the hashed block).

---

## 11. Risks

1. **The deuteron and ⁶Li-config-0 helicity runs refuse at the default
   window** (§3.3). A user who typed nothing but `--isotope d --plan
   helicity-flip` gets a refusal. Mitigation: the message names the cell and
   `--x-max 0.95`; the help and `USAGE.md` list the six. Moving the shipped
   `x_max` is a separate decision (O-2).
2. **The mandatory band's ×2 row needs `--x-max 0.85`** on the default — the
   documented 0.95 recipe breaks for a user who follows it (§3.3).
3. **`--b1-model li6-convolution` alone is refused** under the default (O-1).
4. **The exit-time abort** (§3.4) becomes the default path on any machine
   without `lhapdf.conf` — the pip-wheel case — unless step 1 lands first.
   It is the one step that must not be deferred.
5. **The struck-slot mislabel** (§3.5) becomes reachable from a default
   config; closed by §2 item 6(c).
6. **Caller-supplied kernels** are refused beside the new default under
   today's rule (MEASURED, `test_li7_rank2.py:210-230`); §2 item 6(b).
7. **Unpinned sets** make the default depend on what `/current/` serves on
   the day CI fills its cache (§8 item 5).
8. **The LO caveat reaches every default number** (the development item):
   F₂ᵈ 8.6 % high against NMC, the Bjorken integral without C_Bj, A1 > 1 at
   x = 0.955. Users who read "the shipped default" as "validated" are wrong;
   the banner and the README must say so.
9. **The letter's chain is not the D7 default** (§1.3): numbers quoted from
   LiPolGen at its default are neither the letter's toy footing nor its grid
   column (no EPPS21, no R1998, continuation instead of freeze). Say which
   backend a quoted number used — the rule the tree already states
   ("State BOTH settings next to any plot", `include/lipolgen/pipeline.hpp:619-625`).
10. **Thread safety** of a shared LHAPDF object is tested (§3.6), not proven.
11. **Cost**: ≈ +0.03 s loading and +0.15–0.18 s set-up per run, ×2 per T0/T1
    event; the pytest suite's wall time grows with every `make_config` unless
    the cache (§2 item 5) is adopted — not measured.
12. **Citation-gate churn**: every edited header and `.cpp` moves lines that
    `PHYSICS_CHANNELS.md` cites; `--fix` handles MOVED blocks, EDITED ones need
    a human (step 6).
13. **HepMC3 files do not record the structure-function backend** (only the
    npz `meta` does; `src/hepmc/hepmc_writer.cpp:68-92` writes tool info and
    weight names): a file made before and after D7 differs ×0.80 in
    `xsec_pb` with nothing in it saying why. A run-info attribute is a
    candidate (it touches `docs/HEPMC3_CONVENTION.md`, whose sha256 the T6
    row records) — not part of D7.

---

## 12. Questions only the author can answer (not D7's to decide)

* **O-1 — `--b1-unpol`'s default.** It stays `toy` (`include/lipolgen/pipeline.hpp:955` `B1UnpolSource::Toy`,
  `python/lipolgen/cli.py:821` `b1_unpol="toy"`), and under the D7 default `--b1-model
  li6-convolution` alone is refused by name (§3.3). Keep that (recommended:
  the A = 2 gate's verdict is about `mstw`, so naming the b₁ backend keeps the
  gate's configuration explicit), or let `toy` there follow `--unpol-sf`
  (`li6-convolution` would then fold against CT18NLO, G3b 0.719, "the retired
  stand-in", `python/lipolgen/cli.py:880-882`)?
* **O-2 — the shipped `x_max`.** 1.0 today (`Scenario::x_max`); six
  helicity combinations and the band's ×2 row refuse at the default under D7
  (§3.3). Keep 1.0 and document (recommended for D7), or move the window
  (0.95 cures the six; the band ×2 row needs ≤ 0.85) — a decision of its own.
* **O-3 — the polarized grid's floor.** NNPDFpol11_100 starts at Q² = 1; no
  shipped cell centre is below it (MEASURED), and nothing records it if a
  user's window goes there. Add `LhapdfG1::q2_min()` and
  `meta["pol_sf_grid_q2_min"]` / `["pol_sf_below_grid_frac"]` (the unpolarized
  pair's arrangement) in I7, or leave it?
* **O-4 — the nuclear modification and R.** D7 ships free-nucleon CT18NLO and
  the toy R; the letter's grid chain uses EPPS21 and R1998 (§1.3). An
  `--unpol-sf epps21` option (6Li only; frozen or continued below Q₀) and the
  claim-1 R question (registry row 3; `DECISIONS_PRICED.md` §1) stay separate.
* **O-5 — the letter's wording.** `DECISIONS.md` D7 says "the grid backends
  the letter's production uses"; the letter's published footing is the toy
  (`paper/main.tex:28-31`). §13 proposes wording that says both.

---

## 13. The D7 entry and the development item (texts for `DECISIONS.md`)

The pattern of the file (`docs/open_items/run_2026-10-27/DECISIONS.md:33-48` for D4, `:50-56` for D5/D6;
`../run_2026-09-27/DECISIONS.md` for evidence-bearing entries): a bold
decision sentence, then an italic paragraph with the evidence, what moves and
what does not, and the implementation pointer.

> **D7 — the shipped structure-function defaults become the LHAPDF grid
> backends: `--unpol-sf ct18nlo` (`LhapdfSF("CT18NLO", 0)`) and `--pol-sf
> nnpdfpol` (`LhapdfG1("NNPDFpol11_100", 0)`); `ToyF2` and `ToyG1` stay as
> opt-in options (`--unpol-sf toy`, `--pol-sf toy`); the leading-order
> evaluation of the two NLO sets is filed as the development item below.**
> *(decision 4, on the supervisor's recommendation. The evidence: no backend
> passes its absolute data — NMC F₂ᵈ: toy 8307.9/155, CT18NLO 3401.8/146;
> world g₁ᵈ: ToyG1 107.8/15, 69.4/15, 54.9/28 against NNPDFpol1.1 45.0/15,
> 52.6/15, 51.0/28 (COMPASS, HERMES, E143) — but the grids are fits to the
> world's data where the toys are shapes anchored by eye, CT18NLO passes NMC
> F₂ᵈ/F₂ᵖ (215.5/189, p = 0.090, against the toy's 2161.4/211), and
> NNPDFpol1.1 removes the toy g₁ⁿ's wrong sign over 0.25 < x < 0.6. These are
> the backends of the letter's grid chain (`--pdf grid`, plans/07 WP1), which
> the letter runs and quotes as a bound in its §6 while its published footing
> stays the toy (`paper/main.tex:28-31`); LiPolGen's defaults are that chain
> without its nuclear set (EPPS21), without R1998 and with LHAPDF's
> continuation below Q² = 1.677 instead of a freeze. What moves: every rate
> (⁶Li inclusive σ ×0.79847, coherent ×0.70584), the run-level A_zz ×1.25239
> (exactly the inverse of the rate), inclusive ⁷Li A_∥ ×0.66783; the cos 2φ
> amplitude per cell does not move. The report's `t1_nmc_f2d_over_f2p` row
> flips FAIL → PASS, because its declared rule gives the verdict to the
> shipped default; no other status moves. At the shipped window six deuteron
> and ⁶Li-config-0 helicity runs and the ×2 row of the mandatory b₁ band
> refuse with the cell and the cure (`--x-max 0.95`; 0.85 for the band), and
> `--b1-model li6-convolution` needs `--b1-unpol` named. A missing set refuses
> by name, never falls back. What does not move: the toy physics (pinned with
> the toy named), the polligen port gate and every reference JSON (class-level
> toy objects), no constant, tolerance or registry label. The implementation
> follows `plan_D7_backend_defaults.md`.)*

> **Development item (D7) — the grid sets at next-to-leading order.**
> `LhapdfSF` evaluates CT18NLO, an NLO fit, with the leading-order formula
> F₂ = x Σ e_q² (q + q̄) over five massless flavours
> (`src/lhapdf/lhapdf_sf.cpp:41-46`, `:59-69` `f2_from_weights`), and `LhapdfG1` evaluates
> NNPDFpol1.1, an NLO fit, as g₁ = ½ Σ e_q² (Δq + Δq̄) over d, u, s
> (`src/lhapdf/lhapdf_sf.cpp:47-48`, `:71-82` `g1_scalar`); R stays `r_sigma_lt`. What that costs, measured: NMC
> F₂ᵈ 3401.8/146, 8.6 % high on average, of which massless charm and bottom
> are 3.56 % of the model F₂ᵈ (7.0–7.6 % at x < 0.05) — about half; without
> them 4.6 % high and x-dependent, the rest unexplained, the NLO coefficient
> functions one candidate; world g₁ᵈ fails on all three in-verdict sets; the
> Bjorken integral sits at 98.0–98.2 % of g_A/6 — the C_Bj the LO formula does
> not carry (+19.3 / +12.4 / +9.6 % from g_A/6·C_Bj at Q² = 2 / 5 / 10); and
> the pair's A₁ = g₁/F₁ exceeds the positivity bound at the window's top cell
> (2.02 on the deuteron, 4.41 on the free neutron, x = 0.955), which refuses
> the neutron-tagged and the ⁶Li-config-0 tagged helicity runs (the
> inclusive deuteron's refusal is mostly the tensor term over CT18NLO's
> small F₁ there). The item: evaluate F₂, F_L (hence R)
> and g₁ at NLO in each set's own scheme — the MS-bar coefficient functions
> C_q, C_g and ΔC_q, ΔC_g with Δg, a heavy-flavour treatment consistent with
> CT18 (S-ACOT-χ or FONLL) instead of massless c and b, target-mass corrections
> at low Q² — through an established library (APFEL++, YADISM) or an in-tree
> convolution, and re-run `t1_nmc_f2d`, `t1_hermes_f2d`, `t1_g1d_world` and
> `t5_sum_rules[bjorken]` under their unchanged declared rules. It moves every
> grid-backend number again (σ, A_zz, A_∥), so it is priced before it is
> adopted. Not part of it: the nuclear modification (EPPS21) and R1998, which
> the letter's chain carries and the defaults do not (separate options), and
> the deuteron's own nuclear effects in the F₂ᵈ comparison
> (`validation/benchmarks/t1_nmc_f2d.py:109-119`). Cost not measured.

---

## 14. What this design did not verify

* The C++ suite under design A (the struct cannot be flipped without editing
  `pipeline.hpp`); its counts in §4.4 are from reading, not running.
* `t5_weighted_unweighted` on the D7 kernel; the `--x-max 0.95` cure for the
  four config-0 / config-2 refusals of §3.3 (only the two config-1 deuteron
  cases were re-run); `validation/vmc_tag_fractions.py` on the default.
* The cold-cache set-loading time; the pytest wall time under D7 with and
  without the object cache.
* Anything on GitHub Actions (never run, registry row 14); the wheel step's
  behaviour on the runner.
* Thread safety of a shared LHAPDF object beyond the 1-vs-4-thread
  bit-identity of §3.6.
* Whether D4–D6, implemented first, change any number this plan quotes as
  "today's"; §3's numbers are on `ca9b23a`.
* The sibling's letter text was read, not run; its N_DIS ×0.766 and A =
  −0.2404 are quoted from `plans/07:300-357`. That `moment_A` re-solved on the
  D7 F₁ gives −0.24038 (§3.8) on different inputs (CT18NLO and the toy R at
  ⟨Q²⟩ 4.474, against the letter's EPPS21 and R1998 at 5.069) and is not
  read as a cross-check.

## Appendix — scratch scripts (outside the repository, not permanent)

All under `/tmp/claude-1000/-home-cpeng-Projects-polli-PoliGen/7790c59b-ea07-4677-a36f-0183e2410ef4/scratchpad/d7/`,
run from there after `cd /home/cpeng/Projects/polli/PoliGen && source env.sh`:
`load_time.py` (§3.1 loading), `gen_speed.py` + `gen_speed.jsonl` (§3.1
timing), `confirm_price.py` (§3.2), `matrix.py` (§3.3), `refusals.py`
(§3.3), `delta_bound_d7.py` (§3.8), `run_sim.sh` + `d7pkg/` + `results/`
(§3.7), `harness_*.{base,d7}.txt` (§3.8), `empty_lhapdf/`, `conf_only_lhapdf/`
(§3.4).

## Adversarial check (2026-10-27)

The adversarial check of this design, run on 2026-10-27 against `ca9b23a`
with nothing in the repository edited, reproduced in the checker's own words
(nothing trimmed, no number changed). Its file and line references are the
checker's, at `ca9b23a`, and are written as plain text, not as citations the
records gate checks.
(The `Anchors as of` markers the checker's sixth risk relies on were taken
out of this record on 2026-10-27, after the check, for the reason the
record's opening paragraph gives; its citations are now checked against
the live tree, so the line shifts that risk measures reach this record too.)

### Verdict

The design is sound in its core, but it should not be implemented as written. The recommendation (design A: flip the C++ PipelineConfig selectors and let the LHAPDF tier attach the objects) holds up.

**What reproduces (MEASURED 2026-10-27, ca9b23a, 2000 events, seed 20261026):**
- 6Li sigma 591846.1688784091 -> 472571.91791233147 pb, and 36.18 % of the rate below the CT18NLO grid.
- 7Li A_par for the D7 pair 3.534091e-4.
- Coherent 8806.096 pb and tagged d-p 456479.48 pb.
- The six refusals at x = 0.955.
- The band x2 refusal at x_max 0.95 and 0.9, and the build at 0.85.
- The struck-slot mislabel.
- The abort (exit 134 with no lhapdf.conf) and the clean exit-1 path (lhapdf.conf present, set absent).
- The CI set list and the three tarball sha256s.
- The hashed t4_djangoh block (lines 41-187): the planned edits all lie outside it.

**The record's own NOT-VERIFIED items, now checked:**
- All six refusals build at --x-max 0.95.
- t5_weighted_unweighted PASSES on the D7 kernel (pull +0.359).
- The 6Li alpha RP tag fraction and <k> are bit-identical between toy and D7.
- The F2d/F2p FAIL->PASS also holds on the full Q2 >= 1 cut (245.19/211, p = 0.053).

**Seven defects that block a clean implementation:**
1. **The default depends on the environment.** The D7 numbers depend on whichever lhapdf.conf LHAPDF finds first. Neither set's .info pins Extrapolator, ForcePositive or Interpolator. Changing only lhapdf.conf:
   - Extrapolator nearest: 6Li sigma +5.03 %, coherent +7.02 %, tagged d-p +5.09 %, 7Li A_par x0.952.
   - ForcePositive 1: 7Li A_par x1.488.
   - Extrapolator error: the run is refused.
   
   The pre-flight checks only that the file exists, and meta records none of these settings.
2. **The pre-flight message breaks a string contract.** Its new wording drops LHAPDF's 'Info file not found for PDF set'. _report.py, six harnesses and price_defaults.py use that phrase to label a row 'blocked (environment)'; without it those rows report an error instead.
3. **Two of the design's own rules contradict each other.** The new kernel rule (§2 item 6b) refuses the toy next to a caller kernel. The toy pin plan therefore breaks every caller-kernel test it routes to the toy.
4. **The price table aborts.** The li6-incl-r family of validation/price_defaults.py is refused under D7, which stops the whole regeneration step.
5. **A pinned test file is missing from the pin plan.** test_bench_djangoh_eltail.py pins row 7's verdict and history numbers. Under the recommended 'follow the shipped Born' route the record also plans to move a test bound (0.4 -> 0.474). That conflicts with its own pin rule 2 and the run's 'no tolerance moved' rule; the record's alternative (toy at line 590) should be the recommendation.
6. **Older run records will break under the citation gate.** Simulated D7 line shifts break 89 citations in 19 older run records, yet the record says those records are not edited.
7. **About 20 sites are missing from the inventory.** Among them: the documented C++ examples in USAGE.md, a runtime message carrying the old 7Li sigma, a sha256-pinned data file, and the knob-provenance reason strings, which become false under rule 6b.

### Sites the design missed (21)

| file | line | what |
|---|---|---|
| `validation/benchmarks/_report.py` | 48-51, 195-197 | `LHAPDF_MISSING = "Info file not found for PDF set"` decides blocked (environment) versus error. The pre-flight's new message (record §8 item 1) lacks the substring, so a missing set becomes an error unless the pre-flight keeps LHAPDF's text verbatim. |
| validation/benchmarks/t1_nmc_f2d.py, t1_nmc_f2d_over_f2p.py, t1_hermes_f2d.py, t1_g1d_world.py, t1_r_world.py, t5_sum_rules.py | 366; 375; 375; 453; 464; 504 | `if "Info file not found for PDF set" not in str(exc): raise`. Under the pre-flight a missing set raises (broken harness, exit 2) instead of the blocked sub-row. Same contract in validation/price_defaults.py:155-156 (_ENV_WORDS) and in python/tests/test_bench_t1_nucleon.py:646, test_bench_t1_r_world.py:89, test_bench_nuclear.py:52 and test_bench_report.py:132/:245. |
| `python/tests/test_bench_djangoh_eltail.py` | 40-49, 270-305 (and opt-in 397-425) | Not bench-gated. R_POLRAD_FULL/R_T_PEAK at rel 2e-7, the t-peak shortfall 0.3882 at abs 1e-4 (D7 Born: 1-0.6125 = 0.3875) and HIST_T_PEAK/HIST_POLRAD_FULL ('the 2026-09-23 history numbers ... which the re-wired harness must still reproduce') at rel 2e-7. All fail under the recommended follow-the-Born route. Absent from §5.2 and from the §3.7 simulation. |
| `validation/price_defaults.py` | 563-640 (also 569, 577) | Family li6-incl-r runs b1_model='li6-convolution' with the default b1_unpol. Under D7, validate() refuses it (MEASURED: 'unpol_sf = ct18nlo with b1_unpol = toy on b1_model = li6-convolution'). `_is_environment` is False for 'PipelineConfig...' messages, so line 590 raises and price_table() aborts. Its band rows at X_MAX_CONVOLUTION 0.95 hit the x2 refusal. Line 577 `lg.PipelineConfig().scenario` would trigger the best-effort attach. |
| `src/core/pipeline.cpp` | 2650-2652, 2670-2673 (and 995-1002; python/lipolgen/cli.py:1553-1555) | Knob-provenance reason strings 'which validate() therefore refuses beside it' and 'so validate() refuses any non-toy selector beside one' (plus the clause (3)/(5) comments and the cli comment 'validate() refuses a non-toy selector there anyway'). These become false under §2 item 6b, which accepts ct18nlo beside a kernel and refuses toy. The record marks 2665-2690 KEEP. |
| `python/lipolgen/__init__.py` | 749-751 | Runtime ValueError text: 'measured at the default 7Li beam config 1 and scenario: 590952.42641509 pb for BOTH categories'. Under D7 the default is 470360.06 pb. Pinned by python/tests/test_li7_rank2.py:308 (`assert "590952.42641509" in msg`). Same number at include/lipolgen/pipeline.hpp:1083 and docs/OPEN_ITEMS_SOLUTIONS.md:2792. |
| `docs/USAGE.md` | 140-160, 224-235, 1279-1290, 1713-1725 | The C++ API examples `PipelineConfig cfg; ... Pipeline p(cfg, ...)`. Under design A they throw at validate() (named backend, empty slot) without `#include "lipolgen/lhapdf_sf.hpp"` and `attach_shipped_sf_backends(cfg)`. The li6-convolution example (224-235) also needs b1_unpol named. |
| `docs/USAGE.md` | 2652-2687 | The 'verbatim' KNOB PROVENANCE banner of `lipolgen-run --events 200 --seed 5 --pol-sf nnpdfpol` lists `pol_sf = nnpdfpol` under 'set away from the default'. False under D7. |
| `docs/USAGE.md` | 282-300 | 'the F1 in the denominator is InclusiveKernel's, built from the shared ToyF2 with its R' (Which R each object uses). Under D7 it is CT18NLO's. |
| `validation/benchmarks/data/delta_moment_bound_sources.json` | 164 | "f1": "the shipped 6Li kernel's per-nucleon F1 (ToyF2): the letter's toy column". The file is sha256-pinned (t5_delta_moment_bound.py:303-374 FILE/BODY_SHA256; REPORT.json:4413; REPORT.md:407; validation/benchmarks/README.md:237), so a reword re-pins a vendored digest. The harness docstring also needs rewording: :84-86 'moment_A solved on the shipped kernel's F1', :221, :249. |
| validation/benchmarks/t1_g1d_world.py; validation/benchmarks/t5_sum_rules.py; validation/benchmarks/_plots.py | 659; 672; 34, 50 | PLOT_SLOTS put the toy in the 'default' slot (palette comment: 'slot 1 blue -- the shipped default') with label '(shipped)'. The record lists the equivalent lines for t1_nmc_f2d (:487), t1_nmc_f2d_over_f2p (:505) and t1_hermes_f2d (:606), but not these. |
| `validation/benchmarks/t2_bonus_spectator_shape.py` | 169, 635, 653, 796, 934, 1067 | Row 10's recorded sub-row 'F2n/F2p: the shipped ToyF2.f2n_over_f2p(x)'. Under D7 the shipped species ratio is LhapdfSF::f2n_over_f2p at Q2 = 10. Pinned by python/tests/test_bench_bonus.py:320 (`test_toyf2_hook_is_the_shipped_formula`) and :436. |
| `docs/PACKAGING.md` | 307, 326-333 | Wheel table: `lg.run(channel='inclusive', ...)` 'works' with nothing set, '<x> = 0.027023337580082865 -- bit-identical to the unmasked run'. Under D7 the default run is refused without sets and lhapdf.conf. The record lists only :334. |
| `README.md` | 342-344 | 'As of 2026-10-26 ... 13 pass, 10 fail (recorded), 0 blocked'. The record expects 14/9 after I7. |
| `docs/open_items/run_2026-09-03/STATUS.md` | 81 (row 20), 28 | Registry row 20: '§15.5 D2 decides whether the shipped `toy` default is inside the b1 ban lift'. Same class as row 4 (line 65), which the record does annotate. |
| `docs/benchmarking/08_channel_priority.md` | 36-39, 69, 206, 233, 253 | Ranking text says the shipped inputs are ToyF2/ToyG1 ('a shipped default physics input fails: ToyF2...', 'Rides on. ToyG1', 'F1 from ToyF2 and r_sigma_lt', 'g1p from ToyG1 (--pol-sf nnpdfpol opt-in)') and C2 is scored 'on the toy backend unless a grid number exists'. The record lists only :327, :333, :361-364. |
| docs/OPEN_ITEMS_SOLUTIONS.md; validation/b1_li6_table.py | 525, 1079-1087, 1155; 30, 216, 277-286, 380 | 'the shipped default unpolarised backend, the library ToyF2' (0.440). Same wording the record flags in b1_nuclear.hpp:27-41 and cli.py:902-903; the script prints it in its output. |
| `docs/PHYSICS_CHANNELS.md` | 12-13 | 'Almost every physics ingredient beyond the default toy backends is opt-in' (a citation-gated document). The record lists 120, 128, 136-143, 191, 196, 198-199. |
| `python/lipolgen/__init__.py` | 531 | make_config itself calls `_lipolgen.PipelineConfig()`. Under §2 item 4 (best-effort attach in the binding's init), every make_config, including unpol_sf='toy', loads both grids (about 30 ms each, MEASURED) and touches LHAPDF before the strict attach builds a second pair. The toy escape is then no longer LHAPDF-free. Not addressed. |
| `python/tests/test_knob_provenance.py` | 304-306 (judged at 590-604) | Cell 'b1_band_scale 1 -> 2 on cdks (x_max .95)'. The variant builds on the toy but is refused by positivity under D7 (MEASURED, x 0.871, 1+w_avg -0.008234). It then passes through the `got == "refused"` branch without measuring reach: a cell that keeps passing while testing nothing. |
| `python/tests/test_li6_unpol_band.py` | 5, 31 | Docstrings 'tables on the DEFAULT ToyF2 unpolarised input' and 'shipped default and resolves to ToyF2'. The tests pass under D7 because they name both backends; the text becomes false. |

### Statements found wrong (9)

1. §4.1 row 'src/core/pipeline.cpp:2117, :2665-2690 | knob_provenance ... | KEEP — becomes right by itself once the struct flips'. Only at_default flips by itself. The reason strings at :2650-2652 ('it WINS over --unpol-sf / --pol-sf / --b1-model, which validate() therefore refuses beside it') and :2671-2673 ('validate() refuses any non-toy selector beside one') contradict the record's own §2 item 6b, under which the ct18nlo/nnpdfpol default is accepted beside a kernel and toy is refused.
2. §2 item 6b ('a kernel is accepted beside selectors at their shipped default (compared with PipelineConfig{}) ... a non-default selector beside a kernel is refused as now') contradicts §5.2/§4.4's TOY prescriptions. After D7, Toy is non-default, so a kernel with toy selectors is refused. Sites that break: test_b1_model.py:171-178, :264-267, :795-803 and :940-943, all on `_cfg()`, which §5.2 says to make TOY; tests/test_b1_nuclear.cpp:1988-1990 (`c.kernel = default_inclusive_kernel(LI6()); CHECK_NOTHROW(c.validate()); // ... with the DEFAULT model`) on T11's base(), which §4.4 makes toy_config(). §5.2's 'test_li7_rank2.py:210-230 KEEP ... else TOY' alternative is impossible. The toy_config() replacement is not 'mechanical and greppable' at caller-kernel sites.
3. §5.2 test_bench_chain_rc.py:107-175 RETAKE (recommended): 'the "(v − 1)/(e − 1) < 0.4" bound restated as measured (0.474 ...)'. This moves a test bound after seeing the new value. The re-taken numbers (t-peak 1.01617 -> 1.02051) are the 2026-09-23 history, a dated record, which §5.1 rule 2 says to pin with the toy named. It also contradicts §13's 'What does not move: ... no constant, tolerance'. The alternative (unpol_sf='toy', pol_sf='toy' at t4_djangoh_rad_noRad.py:590) is the one consistent with the record's own rules.
4. §6 and §13, 'Every verdict switch is the harness following its own declared rule, not a new rule'. The declared rules define the verdict configuration operationally: t1_nmc_f2d.py:21-23 'ToyF2() -- the SHIPPED default, i.e. exactly what the generator's `default_inclusive_kernel(deuteron()).tables(x, Q2).f2` returns'; t1_hermes_f2d.py:19-21 and t1_nmc_f2d_over_f2p.py:33-35 say the same. Design A keeps that one-argument kernel the toy, so the rule splits in two. Picking the CT18NLO half, after the 2026-10-26 preview recorded it passing (t1_nmc_f2d_over_f2p.py:86-90), is a rule choice that belongs on the author's question list (O-*). (MEASURED: the flip itself is robust; on the plain Q2 >= 1 cut CT18NLO gives 245.19/211, p = 0.053.)
5. §3.7 '71 failures and 17 errors in 9 files' (and the summary's '9 of 22 files'). The record's own table and scratchpad/d7/results_summary.txt give 11 failing files: test_cli, test_knob_provenance, test_sf_backend, test_b1_model, test_li7_rank2, test_pzz_mode, test_rc_low_x_anchor, test_price_defaults, test_pom_set_band, test_hfs and test_bench_chain_rc.
6. §13 development-item text, '(the inclusive deuteron's refusal is mostly the tensor term over CT18NLO's small F1 there)', also INFERRED in §3.3. MEASURED, inclusive deuteron helicity-flip: at config 2, ct18nlo plus toy g1 builds, while the D7 pair is refused (1+w_avg -0.02576, A1 2.033). At configs 0 and 1 the pair's refused cell (Q2 283.7 / 1007, A1 2.03) is not the tensor-only one (Q2 594 / 1707, A1 0.8446). The vector term (A1 about 2) drives the pair's refusal. The same text's 'neutron-tagged' channel is tagged-d-p, i.e. proton-tagged with a struck neutron (PHYSICS_CHANNELS.md:51).
7. §11 risks, §10 verification and §13 decision text list the six helicity refusals, the band x2 row and li6-convolution, but omit `--inclusive-b1` on tagged-6Li-alpha at the DEFAULT tensor-thirds plan (config 1). It builds on the toy and is refused on the D7 pair (MEASURED: 1+w_avg -0.1302 at x 0.955, Q2 167.3). §3.3 mentions it; the rest of the record drops it.
8. §8 'The rule. Results must not depend on the environment: the default either runs on CT18NLO / NNPDFpol11_100 or refuses by name.' The design does not meet this. §1.1 treats 'LHAPDF's global Extrapolator: continuation (deps/install/share/LHAPDF/lhapdf.conf)' as a property of the tree. It is a property of whichever lhapdf.conf LHAPDF finds first: CT18NLO.info and NNPDFpol11_100.info carry no Extrapolator, ForcePositive or Interpolator, GridPDF.cc:65-69 reads info().get_entry("Extrapolator"), and the pre-flight checks only that the file exists (numbers under the risks).
9. §3.1 table: make_config on the toy measures 0.06-0.07 ms here (three fresh processes, warm cache), not 0.60-0.65 ms. The D7 cost (30-37 ms, second call 16-19 ms), the Pipeline build (0.054-0.062 vs 0.205-0.214 s) and lipolgen-run with 5000 events (0.17 vs 0.37-0.39 s) reproduce.

### Underrated risks (10)

1. **A default that depends on the environment.** MEASURED by changing only lhapdf.conf, same set files, via LHAPDF_DATA_PATH=<dir>::.
   - Extrapolator: nearest: 6Li inclusive sigma 472571.92 -> 496334.16 pb (+5.03 %); coherent 8806.10 -> 9424.48 pb (+7.02 %); tagged d-p 456479.48 -> 479692.86 pb (+5.09 %); 7Li A_par 3.534091e-4 -> 3.364543e-4 (x0.952); below-grid fraction 0.3618 -> 0.3923.
   - ForcePositive: 1: 7Li A_par -> 5.257583e-4 (x1.488).
   - Extrapolator: error: the run is refused ('Point x=0.000288403, Q2=1.0542 is outside the PDF grid boundaries').

   The absolute tensor piece (sigma0 - sigma+-, 276.6234 pb) is the same on all of them, so the published A_zz price is x1.25239 under continuation but x1.19243 under nearest. A pip user must supply their own lhapdf.conf (record §8 item 7), so 36-45 % of every default rate would be set by whoever wrote that file, and meta records nothing about it. Proposed fix: set the extrapolator and ForcePositive explicitly on each PDF in the tier (GridPDF::setExtrapolator, PDF::setForcePositive) and record those settings and the lhapdf.conf path in meta. The +5 % is also the unquantified size of §1.3's 'continuation instead of freeze' difference from the letter's chain.
2. **Pre-flight wording.** If the pre-flight words its refusal as §8 item 1 proposes, a machine without the sets changes report rows from 'blocked (environment)' to 'error' (exit 2), and price_defaults stops recognising the environment (missed sites above). The pre-flight must keep LHAPDF's 'Info file not found for PDF set '<name>'' text. Separately, `lg.lhapdf_quiet()` is bound (python/bindings.cpp:1568) and called directly by t5_delta_moment_bound.py:491, validation/b1_li6_table.py:220 and test_li6_unpol_band.py:74 before any constructor runs. Without lhapdf.conf this still reaches Config::get() and the exit-time abort, whatever the constructors check.
3. **The shipped kernel loses its 1e-12 pin, silently.** T9 (tests/test_b1_nuclear.cpp:1851-1856, 'the only guard on the default b1 path') and test_b1_model.py:335-360 ('pins the SHIPPED default') keep passing on the one-argument toy kernel. After D7 no reference JSON pins the kernel the default run actually builds (CT18NLO F1, so a different Miller-path b1 and Delta), and README:39-41's guarantee is reworded into nothing. A D7 reference behind an LHAPDF guard, or an explicit statement of the loss, is needed.
4. **Tests that pass for the wrong reason:**
   - The knob-provenance cdks x2 cell, which passes through the refused branch.
   - G1b 'default equals ct18nlo/nnpdfpol named', which is true by construction.
   - A G5-type check if cli.DEFAULTS is read from C++; it is tautological unless the literal 'ct18nlo' is asserted.
   - The wheel smoke step, which the record flags. Its proposed replacement, a toy run 'with nothing set', also goes through the best-effort attach. On the runner the vendored libLHAPDF's compiled-in prefix exists, so that 'no data' toy test still loads CT18NLO.
5. **Row 7's verdict depends on the Born.** R moves 1.0199 -> 1.0213 (record §3.8). The header's 'needs no Born at all' (t4_djangoh_rad_noRad.py:29-30) is therefore false at about 1e-3. The recommended route re-takes a pre-registered row's verdict numbers and its 2026-09-23 history after the fact, and the tests pinning them (test_bench_djangoh_eltail.py) are not in the plan.
6. **Citation-gate churn in older records (MEASURED).** I copied the tree to scratch and inserted 4 lines at the first line each D7 edit touches in 18 files. `--records-only` then reports 89 broken citations in 19 dated records, against 0 at baseline. The worst are run_2026-09-03/phase_D_sf_injection.md (21), AUTHOR_DECISIONS.md (11), plan_D4_jm_convention.md (10) and plan_D5_polradfull_default.md (4). plan_D7 itself is safe: its 'Anchors as of' markers exempt every citation in each section (check_physics_channels_links.py:2619-2623). Gate G in §9 steps 6 and 8 cannot pass without editing dated records, which §4 says are not edited. The strict gate's 20 MOVED blocks are handled by --fix, as the record says.
7. **The Python best-effort attach makes the toy path heavy and adds a trap.** It turns every PipelineConfig()/make_config(), including the toy, into a grid load: about 30 ms each, and a second pair of objects on the strict path. If the cli DEFAULTS are read through PipelineConfig() it also happens when the CLI module is imported. Assigning the enum alone (`cfg.unpol_sf = lg.UnpolSfSource.Toy`, a def_readwrite at bindings.cpp:4612) on an attached config is then refused ('toy with an object') — a new trap for users.
8. **Shapes move, not just rates (MEASURED).** <x> of the default tagged-6Li-alpha sample moves 0.02687 -> 0.03408 (+27 %), while the RP tag fraction (0.02372) and <k> are bit-identical. Every default distribution, figure or acceptance-weighted number other than the tag fraction moves, beyond the rate factors the record prices. The upside: the C++ port gate (tests/test_pipeline.cpp:398-472) is SF-independent and needs no toy_config().
9. **The positivity refusal text misleads under D7.** src/core/sampler.cpp:251-264 always says 'e.g. --x-max 0.95 on the shipped grid', but the band x2 row needs 0.85. It also prints 'A1 = g1/F1' on tensor-thirds plans, where g1 is not read (A1 = 0.28 in the band refusal). Under D7 these tensor-driven refusals become common at the default, so the message needs rewording.
10. **The C++ API contract changes.** A default-constructed PipelineConfig no longer runs in a core-only or LHAPDF-off build. The four C++ examples stop being built when the tier is off (CMakeLists.txt:190-200), every documented C++ example refuses, and the proposed conftest stop (exit 3 when the tier is off) makes the LHAPDF tier mandatory for the Python suite. That breaks the project's 'optional tier: suites skip, not fail' convention (ci.yml:255-263) without being stated as a decision.

### The checker's notes: what it ran and confirmed, and what it did not

Nothing in the repository was written. My one `python3 -m lipolgen` timing run left a gitignored build/python/lipolgen/__pycache__/__main__.cpython-310.pyc (17:38); I removed it. `git status` shows only the pre-existing untracked docs/open_items/run_2026-10-27/. Harness imports ran with PYTHONDONTWRITEBYTECODE=1. All scratch work is in /tmp/claude-1000/-home-cpeng-Projects-polli-PoliGen/7790c59b-ea07-4677-a36f-0183e2410ef4/scratchpad/d7adv/: scripts/, the LHAPDF environments lh_{continuation,nearest,error,fp1,empty,confonly}, the gate copies tree/ and tree_shift/, and records_base.txt / records_shift.txt.

**Confirmed by reading:**
- What the backends compute. `--unpol-sf ct18nlo` is free-nucleon Z F2p + N F2n, LO, five massless flavours (lhapdf_sf.cpp:41-69). emc_ratio is set by no pipeline path (pipeline.cpp:1183-1188, 283-288; sf.cpp:127-131). F1 uses r_sigma_lt = 0.18/(1+Q2/50) (sf.cpp:52). nnpdfpol is LO d u s (lhapdf_sf.cpp:71-82). Mixed with toy g1, ToyG1 sits on the kernel's own CT18NLO F1 (xsec.cpp:80-82). pol_sf is not read on tensor-thirds (pipeline.cpp:1338-1351).
- The letter's `--pdf grid` chain: plans/07:263-267, 284-298, 352-357, 1036-1042; paper/main.tex:28-31 and 582-594; structure.py:418, 480-489, 502, 507-508; polligen xsec.py:454-469. All match the record's §1.3.
- CI: ci.yml:73, 79-80, 181-190, 208-223, 240-252, 332-339, 423-438 and build_deps.sh:62, 137-153 match. The tarball sha256s (c9127231…, 48b586b1…, 1afe5f6f…) and member files (375db856…, c64ba5fe…) match the record and the vendored identities in delta_moment_bound_sources.json.
- LHAPDF internals: findFile and paths() never touch Config (Paths.cc:19-53), so the pre-flight idea is sound. The abort mechanism is as described (Config.cc:13-31, Config.h:56-58).
- The hashed DJANGOH block is lines 41-187 (markers at 369-370). Nothing in it names LiPolGen's Born.
- Reference-JSON consumers build class-level kernels (tests/test_reference.cpp:66-99; python/tests/test_kernel_reference.py:29-54), so they are unaffected by design A.
- knob_provenance takes at_default from `const PipelineConfig d;` (pipeline.cpp:2117, 2678, 2690).

**Confirmed by measurement:**
- 6Li sigma, the below-grid fraction and 7Li A_par as stated in the verdict. On tensor-thirds, ct18nlo+toy and ct18nlo+nnpdfpol give identical sigma.
- An independent 66-combination matrix gives exactly the six refusals; all six build at x_max 0.95.
- Band x2: refused at 0.95 and 0.9, builds at 0.85. li6-convolution with b1_unpol toy is refused.
- Missing-set behaviour: no lhapdf.conf -> exit 134 even when Python catches the error; lhapdf.conf present, set absent -> clean exit; the toy run exits 0 in both.
- The struck-slot mislabel reproduces: meta says ct18nlo, sigma = 591846.1688784088.
- Costs and the record-gate simulation as listed under wrong and risks.

**Not verified:**
- The C++ suite under design A (needs pipeline.hpp edited).
- D4-D6 interplay (for example D6's AV18 deuteron moving the tagged-6Li config-0 A1 = 2.02 refusal).
- Thread safety beyond the record's bit-identity test. LHAPDF 6.5.5 has no interpolation cache, only lazy loading (GridPDF.h:208-214).
- Cold-cache load times.
- Anything on GitHub.
- How often real-world lhapdf.conf files differ from the shipped one. I measured the sensitivity, not how often it occurs.
