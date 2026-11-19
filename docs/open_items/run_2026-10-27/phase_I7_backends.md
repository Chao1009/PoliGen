<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Phase I7 — decision D7 implemented: the grid structure functions are the shipped defaults

Run 2026-10-27, branch `claude/run-2026-10-26`, on top of `c52c164` (I6).
Implemented 2026-10-29 per `plan_D7_backend_defaults.md`, amended by its
adversarial check (`checks_1027/D7.json`, reproduced at the end of the plan).
Nothing committed by this phase; a verifier and a fixer follow. Every number
below was measured on this machine on 2026-10-29 unless it says otherwise;
the scripts are in the session scratchpad (`i7/meas/`, `i7/knob/`, `i7/ci/`),
outside the repository.

## 0. What D7 changed, in one paragraph

`PipelineConfig::unpol_sf` defaults to `UnpolSfSource::Ct18Nlo` and
`pol_sf` to `PolSfSource::NnpdfPol`, named ONCE as `UNPOL_SF_DEFAULT` /
`POL_SF_DEFAULT` in `include/lipolgen/pipeline.hpp` and read — never retyped —
by the struct, the bindings (exported), `lipolgen.make_config`, `cli.DEFAULTS`
and the examples; two `static_assert`s tie the texts that call them
`ct18nlo` / `nnpdfpol` to the constants. The core still links no LHAPDF: the
grid objects are **attached lazily** by the LHAPDF tier's filler when a
`Pipeline` is built (`fill_grid_sf_backends`), so `PipelineConfig{}` and
`make_config()` load nothing; the Python module installs the filler at import,
a C++ program calls `install_lhapdf_sf_filler()` once or
`attach_shipped_sf_backends(cfg)`, and the core-only path is the toy NAMED,
`toy_config()`. Every grid the tier opens is **frozen at its knot edges**
(LHAPDF `nearest`, `logcubic`, ForcePositive 0 set ON the object, plus a clamp
of (x, Q²) in `Grid::xfx`), so no `lhapdf.conf` on any machine moves a default
number. A missing set or `lhapdf.conf` is refused **by name** by a pre-flight
that touches no LHAPDF `Config` (exit 1, never the exit-time abort 134), and
keeps LHAPDF's own "Info file not found for PDF set" text. The class-level
null backends and the one-argument `default_inclusive_kernel(ion)` stay the
toy, so the polligen port gate and every reference JSON are untouched; the
toy NAMED is the pre-D7 run bit for bit (24 runs and the kernel, §3.3).
`--b1-unpol` stays `toy`, x_max stays 1.0, EPPS21 and R1998 are not part of
D7. Headline costs: the ⁶Li inclusive σ ×0.83862 (591846.17 → 496334.16 pb),
run-level A_zz ×1.19243; six `helicity-flip` runs, the b₁ band's ×2 row,
`li6-convolution` without `--b1-unpol` and `--inclusive-b1` on the ⁶Li α tag
are refused at their old settings and each refusal computes its cure; one
benchmark verdict flips by its declared rule (`t1_nmc_f2d_over_f2p` FAIL →
PASS).

## 1. The decision, the development item, and what D7 does not move

The decision, as `DECISIONS.md` records it since the fix stage of
2026-11-03 (re-quoted then: until that day this section quoted the `fbc24f2`
entry as "verbatim" although `DECISIONS.md` had been rewritten at `7fd441c`;
the bold text now says "grid sets … (free nucleon, without its EPPS21
ratio)" where it said "grid backends", honesty finding 9):

> **D7 — the shipped structure-function defaults become the grid sets the
> letter's production uses (free nucleon, without its EPPS21 ratio):
> `--unpol-sf ct18nlo` (CT18NLO through LHAPDF) and `--pol-sf nnpdfpol`
> (NNPDFpol1.1); `ToyF2` and `ToyG1` stay as opt-in options; the
> leading-order evaluation of the NLO sets (no coefficient functions,
> massless charm and bottom) is filed as the development item.**
> *(decision 4, on the supervisor's recommendation; implemented 2026-10-29
> per `plan_D7_backend_defaults.md` and `phase_I7_backends.md`. …)*

The author's words (`DECISIONS.md`): *"Your recommendation sounds the best
approach, go with it."* — the author's only words on D7. Below-grid
behaviour is the SUPERVISOR's choice, not the author's (corrected
2026-11-03; this record said "as the author chose it"): the design
(`plan_D7` §1.3) proposed LHAPDF's continuation, its adversarial check
proposed the freeze, and the supervisor took it — the sets are FROZEN at
their grid edges, "LHAPDF extrapolator nearest, set programmatically on every
PDF object the tier builds, so that no lhapdf.conf on any machine can change
a default number". The letter's `--pdf grid` chain freezes Q² only, at Q₀,
on a nuclear set; LiPolGen clamps x and Q² on a free nucleon.

**The development item, in the words the record proposes for it** (§13 of
the plan, with the check's wrong statement 6 corrected): *both shipped grids
are NLO fits evaluated with the leading-order parton formulas — F₂ = x Σ e_q²
(q + q̄) over five massless flavours, g₁ = ½ Σ e_q² (Δq + Δq̄) over three — with
no coefficient functions and massless c and b; 39 % of the shipped ⁶Li rate
sits below CT18NLO's Q² = 1.677 floor, on the frozen edge. Measured costs of
the LO reading: CT18NLO at LO fails NMC F₂ᵈ (3401.8/146, mean model/data
1.086) while F₂ᵈ/F₂ᵖ passes (215.5/189, p = 0.0904); NNPDFpol1.1 at LO fails
the three g₁ᵈ sets (45.0/15, 52.6/15, 51.0/28) and its Bjorken integral sits
12 % above g_A/6·C_Bj at Q² = 5. The six positivity refusals at the default
window's top cell (x = 0.955) are VECTOR-term refusals on the deuteron and
the ⁶Li-α tag (g₁/F₁ ≈ 2 inclusive, 4.4 on the struck neutron of
`tagged-d-p` — proton-tagged, neutron struck); on the inclusive deuteron the
rank-2 term contributes a comparable share (−0.41 / −0.60 / −0.71 of w_avg at
configs 0 / 1 / 2 against the vector −0.59 / −0.42 / −0.31), and with the
toy g₁ configs 0 and 1 are still refused, at a different cell where the
rank-2 term dominates (§3.5).* **Added by the fix stage (2026-11-03,
physics finding 1):** the refusals are the visible half; at LO the shipped
pair breaks |A₁| ≤ 1 from x ≳ 0.83 on the free neutron and from x ≳ 0.87 on
the deuteron and the proton (A₁ⁿ = g₁ⁿ/F₁ⁿ 1.22 at x = 0.85, 1.455 / 1.415 /
1.399 at x = 0.871 for Q² = 10 / 100 / 1000, 1.91 at 0.90; A₁ᵖ 1.006 at
0.871, Q² = 10; the shipped deuteron kernel's A₁ᵈ 0.96 at 0.85, 1.22 / 1.21 /
1.15 at 0.90, 1.53 at 0.93 — the verifier's numbers, the neutron and proton
ones re-taken by the fix stage). `--x-max 0.95` prevents the positivity
REFUSAL but not the violation: the x = 0.871 cell still generates with
|A₁| > 1 on `tagged-d-p` and the inclusive deuteron, not refused because
1 + w_avg stays positive at P_e P_z = 0.49.

**What D7 does not move** (ground rules of the task, each kept): EPPS21 —
the unpolarised input stays the free nucleon Z F₂ᵖ + N F₂ⁿ, and
`EPPS21nlo_CT18Anlo_Li6` stays an option of the validation harnesses; R —
`r_sigma_lt` stays the default, `--r-source r1998` opt-in; `--b1-unpol` —
its default stays `toy` (under the shipped `ct18nlo` the `li6-convolution`
branch must name it, §3.5); the run's x_max — 1.0; no physics constant,
registry label, reference-JSON number, vendored data file or verdict rule
(§6 quotes the rules that moved headlines). One test tolerance DID move after
its value was seen: I7 moved `t5_delta_moment_bound`'s shipped-row error pin
from 1e−15 to `QUAD_REL_TOL`; F7 (2026-11-03) withdrew that: the grid row's
error is recorded, no pass bound; 1e−15 stays on the toy (§13 item 2, §16).

## 2. The design as implemented (design A, amended)

### 2.1 One place

| what | where |
|---|---|
| the two defaults | `include/lipolgen/pipeline.hpp`: `inline constexpr UnpolSfSource UNPOL_SF_DEFAULT = UnpolSfSource::Ct18Nlo;` beside `UnpolSfSource`, `POL_SF_DEFAULT = PolSfSource::NnpdfPol` beside `PolSfSource`; the two `PipelineConfig` fields read them; two `static_assert`s at the end of the header tie the texts |
| Python | `_lipolgen.UNPOL_SF_DEFAULT` / `POL_SF_DEFAULT` exported at the END of the module (after every cited line, as D4 and D6 did); `make_config` resolves `None` to them; `cli.DEFAULTS` derives its two strings from them through `unpol_sf_name` / `pol_sf_name` (not through `PipelineConfig()`, so importing the CLI loads nothing); the `unpol_sf_obj` / `pol_sf_obj` setters re-select the default on `None` |
| the three LHAPDF settings | `LHAPDF_EXTRAPOLATOR` = `"nearest"`, `LHAPDF_INTERPOLATOR` = `"logcubic"`, `LHAPDF_FORCE_POSITIVE` = 0 in `include/lipolgen/lhapdf_sf.hpp`, exported to Python |
| examples | `generate_tagged` / `generate_coherent`: `install_lhapdf_sf_filler()`; `generate_full`: `attach_shipped_sf_backends(cfg)` and the bridge on `cfg.unpol_sf_obj` through `set_pythia_hadronizer` (it drew its T2 species from `ToyF2` under the default before); `generate_inclusive`: its own kernel on `make_lhapdf_unpol_sf(UNPOL_SF_DEFAULT)` / `make_lhapdf_pol_sf(POL_SF_DEFAULT)`. All four run (300 / 300 / 300 / 50 events: 496334 pb inclusive and tagged, 9424.48 pb coherent; the full chain's PYTHIA 50 / 50 ok) |

### 2.2 Lazy attach

`GridSfFiller` is one process-wide function pointer in the core
(`set_grid_sf_filler`, `grid_sf_filler`, under a mutex); the LHAPDF tier's
`install_lhapdf_sf_filler()` installs `attach_shipped_sf_backends`, the Python
module calls it at import. The `Pipeline` constructor calls
`fill_grid_sf_backends(cfg_)` before `validate()` — a no-op beside a
caller-supplied `kernel` (the selectors are not read there), when no grid
slot is empty, and without a filler. `validate()` accepts an empty
`unpol_sf_obj` under `Ct18Nlo` (or `pol_sf_obj` under `NnpdfPol`) only while a
filler is installed, and otherwise refuses naming the tier, the set, the two
helpers and the toy escape. The objects are ONE per process
(`make_lhapdf_unpol_sf` / `make_lhapdf_pol_sf`, built on first use under a
mutex, warmed before sharing) and shared read-only; threaded generation on
them is bit-identical to single-threaded (new C++ case). **`b1_unpol` is NOT
lazy** (changed during implementation): an empty `Ct18Nlo` b₁ slot is
refused by `validate()` exactly as before D7
(`test_validate_refuses_a_named_backend_with_an_empty_slot` pins that
contract), and `set_b1_unpol(Ct18Nlo)` attaches the same cached object
eagerly — so under the shipped `--unpol-sf ct18nlo` a `--b1-unpol ct18nlo`
run folds against THE SAME object as the rate, and the partial cancellation
of the two F₁ is exact again. `set_unpol_sf` / `set_pol_sf` pre-flight the
grid and leave the slot empty for the `Pipeline`.

### 2.3 The freeze

`Grid::xfx` clamps (x, Q²) to the grid's own knot edges and evaluates;
`Grid::open` also sets `setExtrapolator("nearest")`, the interpolator
(`logcubic`, only when it differs) and `setForcePositive(0)` on every
`GridPDF` the tier builds — `LhapdfSF`, `LhapdfG1` and both grids of
`Epps21Ratio`. The clamp IS LHAPDF's `nearest` to the bit (new C++ case,
against a raw `GridPDF` with `nearest` set): `NearestPointExtrapolator`
re-interpolates at the clamped knot but copies the whole `KnotArray` per
call, 3.7 s against 0.016 s for 5000 events, so the clamp is what runs and
the setter is the documented equivalent. Measured consequences:

* **No `lhapdf.conf` moves a default number.** Two fresh processes — this
  store (its conf: `logcubic`, `continuation`, `ForcePositive: 0`), and the
  same set files under an `lhapdf.conf` saying `continuation`, `linear`,
  `ForcePositive: 1` (`LHAPDF_DATA_PATH=<dir>::`; the SAME extrapolator, so
  this pair tests the interpolator and ForcePositive, and the extrapolator
  rests on the absolute σ pin — named so 2026-11-03, honesty finding 5; the
  fix stage added a variant saying `Extrapolator: error`,
  `test_the_default_numbers_do_not_depend_on_the_conf_extrapolator`) —
  give the same σ and per-category σ to the bit on ⁶Li inclusive, coherent,
  `tagged-d-p` and ⁷Li inclusive (`test_the_default_numbers_do_not_depend_on_lhapdf_conf`,
  new). The pre-D7 build under the same two files (interpolator and
  ForcePositive differing, extrapolator not): ⁷Li σ 470360.0619 →
  470574.0414 pb (+4.5 × 10⁻⁴). A different comparison — this store's conf
  against one saying `Extrapolator: nearest`, i.e. continuation against
  the freeze — moved the pre-freeze ⁶Li σ by +5.03 % (`plan_D7`, underrated
  risk 1); `DECISIONS.md` quotes that one.
* `meta` records it: `unpol_sf_extrapolation` / `pol_sf_extrapolation`
  (`"nearest"`; `"none (...)"` on the toy), `pol_sf_grid_q2_min` (1.0) and
  `pol_sf_below_grid_frac` (0.0 on the shipped window) beside the existing
  `unpol_sf_grid_q2_min` (1.677025) and `unpol_sf_below_grid_frac`; the
  banner prints the floor, the fraction, the word FROZEN and the LO caveat.
  Since F7 (§16, V7-3) the three POLARISED keys and the banner's polarised
  line are written only where g₁ RAN; on a run that does not read it (every
  tensor plan, the CLI default; the coherent channel) the extrapolation key
  carries `meta["pol_sf"]`'s label and the two numbers are NaN — I7 credited
  `nearest`, 1.0 and 0.0 to a backend that never ran.
* The bound side of `t5_delta_moment_bound` reads `Epps21Ratio`, whose grid
  starts at x = 1e−7 while the harness integrates from 1e−9: ⟨x⟩_g^EPPS21
  moves −4.8e−8 / −1.8e−6 / −1.9e−6 relative at Q² = 2 / 5 / 10 (the freeze
  below x = 1e−7), and with it the B_tight edge (`saturation_tight` of
  `moment_A`, re-taken). B, B_loose, ⟨x⟩_g and ⟨x⟩_Δg do not move.

### 2.4 The pre-flight (amendment 2)

`lhapdf_set_problem(set, member)` uses `LHAPDF::findFile`,
`findpdfsetinfopath`, `findpdfmempath` and `paths()` only — no `Config` is
constructed, so the refusal path cannot reach the exit-time abort — and
returns, on a miss: *"LHAPDF cannot open the set 'CT18NLO': Info file not
found for PDF set 'CT18NLO'. LiPolGen's pre-flight (D7, 2026-10-27), run
before any LHAPDF call: searched …; LHAPDF_DATA_PATH=…; lhapdf.conf NOT found
/ found at …. Install the set (…build_deps.sh…, sha256-pinned) … or run the
toy explicitly with --unpol-sf toy (Python unpol_sf='toy', C++
UnpolSfSource::Toy or lipolgen::toy_config()) and say so. It is never
substituted by a toy in silence."* `_report.py`'s `LHAPDF_MISSING`, the six
harnesses' `"Info file not found for PDF set"` contract and
`price_defaults.py`'s `_ENV_WORDS` still recognise a missing set (the
substring is kept) — **but, as I7 left them, not where it mattered** (this
sentence said they "keep working"; corrected by F7, §16, V7-1): since D7 the
four T1 rows' VERDICT configuration and every price family's DEFAULT row
are the configurations that need the sets, and with a set missing the four
T1 harnesses crashed (KeyError `'dis'` / `'experiments'`, exit 2, `error`
in the report, where c52c164 reported the toy's verdict with the grid
sub-rows blocked) and three price families did too (li6-incl, li6-coh and,
with only NNPDFpol11_100 missing, li6-tag). F7: the T1 rows raise LHAPDF's
own words (`blocked (environment)`, exit 2 — the rule of rows 4 and 13) and
the recipe marks every probe of an `environment` default row `environment`.
With the set present and only `lhapdf.conf` missing, the advice now says to
supply `lhapdf.conf` rather than to install the set (V7-10).
`lhapdf_quiet()` now calls `setVerbosity` only beside an `lhapdf.conf`.
Re-enacted: `lipolgen-run --events 10` under an empty data path and under a
data path holding only `lhapdf.conf` exits **1** with that text; the toy
named exits 0; pytest stops with return code **3** naming both sets
(`conftest.py`); a Python `lg.run()` raises the `RuntimeError` and the
interpreter exits 0 after catching it.

### 2.5 The configuration rules (amendment 3)

* **Kernel beside selectors, one rule**: accepted at the SHIPPED values
  (compared with `PipelineConfig{}`; `meta` says "caller-supplied kernel"),
  refused at any other — `toy` included — with a message naming both
  (`src/core/pipeline.cpp`, `validate()` clause (3)). The knob-provenance
  reason strings at src/core/pipeline.cpp lines 2650–2652 and 2671–2673 and the comment
  at python/lipolgen/cli.py lines 1553–1555 were rewritten to say this (they said a non-toy
  selector is refused beside a kernel).
* **Struck slots**: a directly set `struck.f2_source` / `struck.g1_model` is
  refused unless the selector is `Custom` (it was refused under `Toy` only;
  measured before the fix, a `tagged-6Li-alpha` run with `ToyF2` there and
  `unpol_sf = ct18nlo` ran at 591846.1688784088 pb, the toy's, labelled
  ct18nlo).
* **The T2 species draw**: `check_bridge_species_backend` moved into the C++
  PYTHIA tier (`src/pythia/bridge_config.cpp`), so a C++ caller binding a
  default-option bridge to a default config is refused too (the example did
  exactly that); the Python `set_pythia_hadronizer` fills the grid slots
  first and delegates. `lipolgen-run` and `lipolgen.run` attach and hand the
  bridge `cfg.unpol_sf_obj`.
* **The positivity refusal** (`src/core/sampler.cpp`, line-neutral): the
  message now splits w_avg into the VECTOR (P_e P_z D(y) A₁) and RANK-2 parts,
  counts the refused cells of that spin state, names the lowest x refused
  and COMPUTES the `--x-max` cure ("below x = 0.871" for the band's ×2 row,
  not the fixed "--x-max 0.95" it printed).

## 3. Measured

### 3.1 The numbers that move (6Li config 1 unless stated; σ is the cell integral, N-independent)

| quantity | toy (the pre-D7 default) | D7, frozen | the plan's D7 (continued, pre-freeze) |
|---|---|---|---|
| ⁶Li inclusive σ | 591846.1688784091 pb | **496334.155954048** (×0.838620) | 472571.9179 (×0.79847) |
| ⁶Li coherent σ | 12476.025113185518 | **9424.476321023749** (×0.755407) | 8806.096 |
| d `tagged-d-p` σ | — | **479692.8570609205** | 456479.48 |
| ⁷Li inclusive σ (helicity-flip) | 590952.4264150970 | **494062.76099941705** (×0.836045) | 470360.0619 |
| ⁷Li A_∥ (helicity-flip, P_e P_z 0.49) | — | **3.364543e−4** | 3.534091e−4 |
| run-level A_zz, ⁶Li inclusive | −5.193231e−4 | **−6.192590e−4** (×1.19243) | ×1.25239 |
| A_∥, ⁶Li inclusive | −1.171716e−3 | **−1.515385e−4** (×0.12933) | ×0.13583 |
| below-grid fraction of the run's own rate: inclusive / coherent / d–p / ⁷Li | 0 (no grid) | **0.3923453 / 0.4837146 / 0.3951002 / 0.3928570** | 0.3618 / 0.4475 / — / — |
| polarised floor; fraction below | — | Q² = 1; **0** | — |

The frozen numbers are the plan's checker's `nearest` numbers to every
printed digit (its 496334.16, 9424.48, 479692.86, 3.364543e−4, 0.3923), as
they must be: the freeze IS `nearest`.

### 3.2 Environment independence — §2.3.

### 3.3 The toy NAMED is the pre-D7 default, bit for bit

`i7/meas/m_toy_bitident.py` run in the c52c164 build and in this tree, every
run with `unpol_sf="toy", pol_sf="toy"`: σ, the per-category σ and an md5 over
every numeric generated column of 24 runs (6Li inclusive / coherent /
tagged-α, 7Li inclusive / tagged-α, d inclusive / d–p, ³He inclusive, each at
configs 0 / 1 / 2) plus `default_inclusive_kernel(li6()).tables(0.3, 10)` —
**25 of 25 identical**. The C++ D7 case pins `toy_config()`'s 6Li σ at
591846.1688784091 too.

### 3.4 Cost (amendment 4; fresh process each, median of 3, idle machine)

| | c52c164: toy (default) | c52c164: grid named | D7: default | D7: toy named |
|---|---|---|---|---|
| `make_config()` (first / second call) | 0.064 / 0.009 ms | **36.8 / 15.6 ms** (eager grid load) | **0.14 / 0.036 ms** (nothing loaded) | 0.070 / 0.012 ms |
| first `Pipeline` build / second | 0.061 / 0.048 s | 0.206 / 0.208 s | 0.248 / 0.207 s | 0.054 / 0.053 s |
| `generate` 5000 events | 4.4 ms | 11.8 ms | 9.4 ms | 4.5 ms |
| `lipolgen-run --events 5000` wall | 0.17 s | 0.36 s | **0.38 s** | 0.16 s |

`import lipolgen` 0.06–0.07 s in all four (the filler installs nothing). The
D7 default costs the first `Pipeline` ≈ 0.19 s of grid load and table
evaluation, once per process; no grid is loaded by `make_config()`,
`PipelineConfig()` or importing the CLI (amendment 4: `cli.DEFAULTS` is not
read through `PipelineConfig()`).

### 3.5 Where the defaults are refused (amendment 6)

`i7/meas/m_refusals.py`: all 66 isotope × configuration × channel × plan runs
at x_max 1.0, the default against the toy named. The toy refuses none; the
default refuses six, every one `helicity-flip` at P_e 0.7, at x = 0.955, and
every one builds at `--x-max 0.95`:

| run | Q² | 1 + w_avg | vector (A₁) | rank-2 | cells | ct18nlo + toy g₁ | toy F₂ + nnpdfpol |
|---|---|---|---|---|---|---|---|
| d inclusive, config 0 | 283.7 | −0.004833 | −0.5914 (2.026) | −0.4134 | 10 | refused (Q² 594: 1 + w −0.04631, vector −0.5259 at A₁ 0.845, rank-2 −0.5205) | builds |
| d inclusive, config 1 | 1007 | −0.01952 | −0.4186 (2.024) | −0.6009 | 7 | refused (Q² 1707: −0.007475, −0.3135 at 0.845, −0.6939) | builds |
| d inclusive, config 2 | 1897 | −0.02576 | −0.3126 (2.033) | −0.7132 | 1 | builds | builds |
| `tagged-d-p`, config 0 | 229.7 | −0.005726 | −1.006 (4.366) | 0 | 12 | builds | builds |
| `tagged-d-p`, config 1 | 1119 | −0.02414 | −1.024 (4.406) | 0 | 6 | builds | builds |
| `tagged-6Li-alpha`, config 0 | 480.9 | −0.02789 | −1.028 (1.98) | 0 | 5 | builds | builds |

**Not only refusals** (fix stage, 2026-11-03, physics finding 1): the
`--x-max 0.95` cure removes negative densities, not unphysical asymmetries —
the shipped pair has |A₁| > 1 from x ≳ 0.83 (neutron) and x ≳ 0.87
(deuteron, proton), so the x = 0.871 cell of the d inclusive and
`tagged-d-p` rows builds at 0.95 with |A₁| > 1 (§1).

**The inclusive-deuteron split, as measured** (the check's wrong statement 6):
the default pair's refused cell is NOT the tensor-only one; both terms are
negative and of comparable size there; with the toy g₁ configs 0 and 1 still
refuse, at another cell where the rank-2 term dominates; with the toy F₂ none
refuses. The tagged rows are the vector term alone. Three more recipes:

* **the b₁ band's ×2 row** (`--b1-band-scale 2`): `cdks` refused at x_max 0.95
  and 0.9 (x = 0.871, Q² = 186, 1 + w_avg = −0.008234, the rank-2 term alone,
  23 cells), `li6-convolution` on `--b1-unpol ct18nlo` and `mstw` likewise
  (Q² 315.3; −0.005282 / −0.0109); all three build at **0.85**; bands 0 and 1
  build at 0.95; the toy builds the ×2 row at 0.95;
* **`li6-convolution` with `--b1-unpol` left at `toy`** is refused by
  `validate()` under the default ("unpol_sf = ct18nlo with b1_unpol = toy …
  the shipped default since D7; name --b1-unpol: --b1-unpol ct18nlo to fold
  against the same one, --b1-unpol mstw for the configuration the A = 2 gate
  passes on (G3b 0.843243), or --unpol-sf toy for the pre-D7 toy pair"); named
  `ct18nlo` it builds;
* **`--inclusive-b1` on `tagged-6Li-alpha`** at the default tensor-thirds
  plan, config 1: refused (x = 0.955, Q² = 167.3, 1 + w_avg = −0.1302, rank-2
  −1.13, 24 cells); builds at `--x-max 0.95`; the toy builds it at 1.0.

### 3.6 Shapes move, not only rates (the check's risk 8)

20 000 events, seed 20261026, tensor-thirds: ⟨x⟩ 0.02667 → **0.03278**
(inclusive), 0.02674 → **0.03287** (`tagged-6Li-alpha`), 0.00489 → 0.00524
(coherent); ⟨Q²⟩ 4.18 → 4.54 / 4.18 → 4.47 / 2.65 → 2.86 GeV². The RP tag
fraction (0.03155) and ⟨k⟩ (0.12229 GeV) of that tagged sample are identical
to five digits; at 400 k events the tagging-optics fraction moves 0.248635
→ 0.2486425 and ⟨k⟩ in the sixth digit (§5).

### 3.7 Top-cell b₁/F₁ on the shipped kernel

CT18NLO's F₁ at the window's top cell (x = 0.954993, Q² = 167.344) is ×0.134
of `ToyF2`'s, so b₁/F₁ there is **48.67** (`cdks`), **42.85**
(`li6-convolution` on `--b1-unpol ct18nlo`), **1.085** (`miller`) against the
toy's 6.52 / 5.91 / 0.145. `cli.B1_TOP_CELL_B1_OVER_F1` (the one definition
the `--x-max` help and the `--b1-model` refusal read) and its test were
re-taken on the shipped kernel; the toy values are in its comment.

## 4. Every re-taken number, before → after (tests; each with a dated comment naming D7)

**C++.** `tests/sf_test_support.hpp` (new) carries the one rule — `toy_config()`
for every machinery test that builds a `Pipeline` (the core suite installs no
LHAPDF filler) — applied to the root construction sites of
`tests/test_pipeline.cpp`, `test_b1_nuclear.cpp`, `test_rc_pipeline.cpp`
(`base_cfg`), `test_t2.cpp`, `test_breakup.cpp` (`tagged_cfg`), `test_fsi.cpp`
and `test_triton_sf.cpp`, each with the dated comment. Re-taken: the D2 case
(the shipped default is refused without a filler; a kernel beside the shipped
selectors accepted, beside `toy` refused), T11 / T16's two caller-kernel
sites in `test_b1_nuclear.cpp` (selectors set to the shipped values beside
the kernel, so the b₁ guard and not the SF guard is what refuses). New, in
`tests/test_lhapdf.cpp`, 5 cases: the pre-flight in both missing-set
environments (LHAPDF's own words, the toy escape, no `Config`); the freeze
equals LHAPDF's `nearest` to the bit below both floors; the default needs
the filler, with it runs on the ONE shared pair at 496334.155954048 pb, the
toy named at 591846.1688784091 pb; the kernel-beside-selectors and
struck-slot rules; threaded generation bit-identical on the shared grids.

**Python** (the files that changed; every other file untouched):

| file | what | before | after |
|---|---|---|---|
| `test_sf_backend.py` | G1b: the default asserted by the LITERAL names and run against the same names spelled out (no tautology) | default = toy | default = ct18nlo / nnpdfpol; the toy named = the pre-D7 σ |
| … | tagged σ ratio, `meta["unpol_sf_below_grid_frac"]`, banner | 0.7984709; 0.3617907; "36.18 %" | **0.8386202**; **0.3923453**; "39.23 %", FROZEN, the LO caveat |
| … | new meta keys | — | `*_extrapolation` = nearest, `pol_sf_grid_q2_min` 1.0, `pol_sf_below_grid_frac` 0 |
| … | the pre-flight, the kernel and struck-slot rules, `cli.DEFAULTS` literal | — | pinned |
| `test_d7_shipped_kernel.py` (new) | the shipped 6Li kernel at six (x, Q²) points (one below the floor) and the four sweet spots at rtol 1e-12; = the default `Pipeline`'s own kernel; frozen below the floor; the opposite-`lhapdf.conf` bit identity | — | 5 tests |
| `test_b1_model.py` | `_cfg` names the toy pair (the b₁ selectors' tests pin toy physics); `cli.B1_TOP_CELL_B1_OVER_F1` on the shipped kernel | 6.52 / 5.91 | **48.67 / 42.85** (and two D7 tests: `li6-convolution` alone refused; the band ×2 row at 0.85) |
| `test_knob_provenance.py` | SF cells `ct18nlo -> toy`, `nnpdfpol -> toy`; the kernel cell on `ToyF2` (a CT18NLO kernel is bit-identical to the default on coherent); `CONV` names `--b1-unpol ct18nlo`, its b1_unpol cell `ct18nlo -> mstw`; the cdks ×2 cell at x_max 0.85; `inclusive_b1` at 0.95; the T2 cells wire the bridge | — | 627 of 627 cells run with the same status class as on c52c164 (§7) |
| `test_li7_rank2.py` | the selector-reach ratios (toy against the default) and the `make_plan` refusal's measured zero | 0.795936, 0.792543; "590952.42641509" | **0.836045, 0.832599**; "494062.76099942" (the toy's beside it in the text) |
| `test_pzz_mode.py`, `test_rc_low_x_anchor.py` (A2), `test_pom_set_band.py` | sec. A3 / A2 / D4 numbers measured on the toy | — | the toy NAMED; D7 values in the comments (A_zz −5.124064e−4 … at the A2 points) |
| `test_price_defaults.py` | the default row; `RATIOS` | σ 591846.1688784091, A_zz −5.193231e−4, A2_TOY; mstw ×0.7934 / ct18nlo ×0.7985 | σ **496334.155954048**, A_zz **−6.192590e−4**, A2 **−6.520872 / −5.317408 / −7.145333 / −6.633634e−4**; mstw ×0.9460 / **toy ×1.1924**; the toy ROW pins the old numbers (cited lines 38–44 untouched) |
| `test_cli.py` | `--isotope D --plan helicity-flip` | x_max 1.0 | `--x-max 0.95` for every spelling (the deuteron's refusal) |
| `test_hfs.py` | the bridge | own ToyF2 | `fill_grid_sf_backends` + `f2_source = cfg.unpol_sf_obj` |
| `test_bench_t1_nucleon.py` | headlines and the row statuses | toy, FAIL ×3 | ct18nlo / ct18nlo / nnpdfpol; FAIL / **PASS** / FAIL; CT18NLO's `all` χ² 3523.9906 → **3959.3713** (the points below the floor, now frozen) |
| `test_bench_t5_delta_bound.py` | the shipped row (CT18NLO's F₁), positivity, the bound's EPPS21 edge, the REPORT line | M2 1.238347e−4 …, w_min 0.9493997990, A_T [−0.1012, 0.0506] | M2 **1.460361e−4 / 1.481675e−4 / 1.501341e−4**, saturations re-taken, w_min **0.9134944997** (x = 0.9, Q² = 10), A_T **[−0.1730, 0.0865]**; the shipped row's numerical-error pin MOVED after its CT18NLO values were seen — from `e < 1e-15` absolute (which 3.3e−14 / 4.9e−14 / 5.5e−14 fail) to the harness's declared `QUAD_REL_TOL` (~10⁴ looser) — and F7 (2026-11-03) withdrew that bound: the grid row's error is RECORDED with no pass bound, the 1e−15 pin applies, unchanged, to the TOY integrand (`letter_kernel()`) (§13 item 2, §16) |
| `test_bench_t5_closures.py` | the weighted/unweighted closure on the shipped kernel | N 3032577 / 2987706, σ 89.193441 / 89.167357, shift 1.4717e−2, no-op +18.29 | N **2712940 / 2669772**, σ **79.792353 / 79.767651**, err_w 0.0488204, pull +0.359, σ_cat 79.767466, σ_unpol 78.522871, shift **1.58501e−2**, no-op **+18.61** |
| `test_doc_link_gate.py` | the as-of census | 110 | **119** (the nine record citations D7's edits moved, §9) |
| `test_bench_bonus.py`, `test_li6_unpol_band.py` | docstrings ("shipped" ToyF2) | — | labelled |
| `conftest.py` | — | — | the session check (return code 3 without the tier or the sets) |

## 5. The numbers I4, I5 and I6 re-took on the toys, re-measured on D7

I4, I5 and I6 re-measured their sample-level numbers on the structure
functions shipped then — the toys — and each record said D7 would move them
again. Their own scripts (`i6/meas/m_fsi.py`, `m_fsi_pct.py`, `m_rc.py`,
`m_tagfrac.py`, `m_t1_virt.py`; `i5/meas/m_points.py`, `m_extra.py`,
`d5chk/clip_configs.py`), re-run unchanged on the D7 tree except that the
band half-widths the I5 scripts had hard-coded from the toy kernel are now
taken from the shipped kernel by the same recipe (which reproduces the toy's
constants exactly when run on the toy). The toy columns are the I4 / I5 / I6
records' numbers; the toy-named path is the pre-D7 run bit for bit (§3.3), and
where it was re-run here (the tag fractions, the T1 virtuality, the RC tail at
Q² = 2 / 5 / 10 and the whole-run ⟨rc_tail⟩, the band half-widths) it
reproduces them to the last printed digit.

| quantity (recipe) | toy SFs (I4 / I5 / I6) | **D7** | where it is published, labelled |
|---|---|---|---|
| N2 sweet-spot a₂ (four cells) | the I4 values | **identical to 1e−14** (F₂ cancels at fixed Δ/F₁) | `test_d7_shipped_kernel.py` pins both |
| A2(S1…S4) per unit P_zz, `price_defaults` super-bins | −6.514569 / −5.312811 / −7.145394 / −6.645439e−4 | **−6.520872 / −5.317408 / −7.145333 / −6.633634e−4** (through the cells' weights only) | `test_price_defaults.py`, the price table |
| t5 positivity w_min, A_T | 0.9493997990 (x 0.9, Q² 2), [−0.1012, 0.0506] | **0.9134944997** (x 0.9, Q² 10), **[−0.1730, 0.0865]** | `test_bench_t5_delta_bound.py`, `validation/benchmarks/README.md` |
| RC band clamp, the published recipe (⁶Li seed 1 / seeds 1–5; ⁷Li seed 11 / seed 1) | 16 / 16, 20, 21, 30, 34; 507 / 499 | 16 / 16, 20, 21, 30, 34; **508** / 499 (Hulthén: 245 / 245, 223, **268**, 219, 235; 520 / 502) | USAGE §7b (2.540 % on ⁷Li) |
| ⁶Li α D3 FSI, VMC: Σw cluster / nucleon; Σw/Σw_off | 15989.35 / 16428.71; 0.799467 / 0.821436 | **15989.24 / 16428.58; 0.799462 / 0.821429**; survival, w ranges, percentiles, 64.515 %, clipped cells unchanged; Hulthén unchanged | USAGE §3 |
| ⁷Li α and d–p FSI tables | — | unchanged to every printed digit | — |
| ⁶Li α tag fraction (400 k, YR / tagging-legacy / tagging), ⟨k⟩ | 0.0334375 / 0.287445 / 0.248635, 0.12264865 GeV | 0.0334375 / 0.287445 / **0.2486425**, **0.12264912** | USAGE §5 |
| ⁷Li α tag fractions (configs 0 / 1 / 2) | 0.99846 / 0.99846 / 0.99846 (tagging 0.99955 / 0.99938 / 0.99951) | unchanged | — |
| T1 ⟨p_N² − M_N²⟩ (50 k), ⁶Li α / ⁷Li α / d–p; redraws | −0.05817 / −0.08468 / −0.03609; 27 / 7 / 3 | **−0.05811 / −0.08480** / −0.03609; **32** / 7 / 3 (0.064 %) | `docs/T2_CHAIN.md`, USAGE §3 |
| RC tail, default `polrad-full`, x = 0.01 / 0.03 / 0.10 / 0.30, Q² = 5: w_tail − 1 | 6.62870e−4 / 1.75050e−4 / 1.95963e−4 / 4.88356e−4 | **7.50159e−4 / 1.85041e−4 / 1.78306e−4 / 3.58654e−4** | `OPEN_ITEMS_SOLUTIONS.md` §9 |
| … at Q² = 2 / 10 | 1.71661e−3 / 2.19742e−3; 2.52430e−3 / 2.51627e−5; — / 5.04628e−5 | **2.26124e−3 / 2.39666e−3; 2.25917e−3 / 2.34919e−5; — / 4.01528e−5** | same |
| … ΔA_zz of the whole tail | −1.178978e−7 / +9.821438e−8 / +2.995857e−7 / +2.378230e−8 | **−8.874874e−8 / +1.118125e−7 / +2.480364e−7 / +1.282748e−8** | OPEN_ITEMS §9 |
| … the band's lead at x = 0.01 (anchors 0.30 / 0.266 / 0.113; `a_transfer_frac` 0.5 / 1); at x = 0.10 | ×1152.2 / ×1021.7 / ×434.0; ×1288.3 / ×1629.5; ×323.1 | **×1732.1 / ×1535.8 / ×652.4; ×1936.5 / ×2449.6; ×355.1** | PHYSICS_CHANNELS 384–385, OPEN_ITEMS 813, 907, USAGE 2277–2289, `rc.hpp` |
| … the knob spreads at x = 0.01 (qe_suppression, qe_kf_gev, c0_shape, qe_tensor_scale, fq_scale, tail_tensor_scale) | 1.165e−7, 5.819e−8, 3.467e−8, 2.650e−7, 8.165e−7, 3.656e−8 | **1.491e−7, 7.450e−8, 3.894e−8, 2.998e−7, 9.240e−7, 4.137e−8** — the ordering and "`qe_kf_gev` the smaller QE knob" hold | OPEN_ITEMS 836–844 |
| … `qe_tensor_scale` 1 as % of the band at x = 0.01 / 0.10 / 0.30 | 0.195 / 0.018 / 1454 % | 0.195 / 0.018 / 1454 % (both scale with A_zz) | USAGE 2249 |
| … ΔA_zz(vmc-ft)/ΔA_zz(ho) − 1 | +29.4 / −0.0 / +0.0 % (t-peak +18.7 / −4.0 / +90.7) | **+43.9** / −0.0 / +0.0 % (t-peak +22.5 / −4.6 / +109.7) | USAGE 2267 |
| … mean `rc_tail`, 5000-event CLI files (P_z 0.7) | t-peak 1.020341, polrad-full 1.028010 (max 3.2380 / 4.0053) | **1.027211 / 1.037150** (max 3.8526 / 4.8560), 0 clipped — recipe (seed) not kept; re-taken 2026-11-03 with the CLI default seed (`lipolgen-run --events 5000 --rc tensor-band --rc-tail-model t-peak` / `polrad-full`): **1.027943 / 1.037955** (max 3.7617 / 4.6656) | — |
| … whole-run mean `rc_tail`, 200 k, seed 1234, t-peak / t-peak+ll / polrad-full, P_z 0.7 | 1.021778527 / 1.040274188 / 1.029702912 | **1.025494331 / 1.047116364 / 1.034783841** | PHYSICS_CHANNELS 828–831, OPEN_ITEMS 832, USAGE 2449; USAGE 2377–2379, `rc.hpp`'s `tail_model` note and the `--rc-tail-model` help labelled only by the fix stage, 2026-11-03 (code finding 3: this row named them as labelled) |
| … at P_z 0 | 1.021836305 / 1.040368933 / 1.029775347 | **1.025506240 / 1.047148070 / 1.034803052** | same |
| … node clip at `tail_max` on configs 0 / 1 / 2 | 0.3986 / 0.3858 / 0.3600 % | unchanged (model-level) | — |

The ratios of tail pieces (σ^q_U/σ^el_U, the QE shares, (1/6)σ^el_T/σ^el_U,
the Pauli cut QRT(S)/QRT(1)) do not move, as I5 predicted: they carry no
Born. The tagged-6Li-α `--inclusive-b1` A_zz the documents quoted was the
Hulthén-era toy number (−1.557969e−3); re-taken on D7's default (VMC wave,
`--x-max 0.95`): −1.660420e−3 (toy) → **−1.979943e−3** (ct18nlo, ×1.19243) →
−2.092855e−3 (mstw, ×1.26044).

## 6. Benchmark verdicts: the declared rules, quoted, and what moved (amendment 8)

A verdict was re-pointed ONLY where the harness's declared rule, quoted
here from `c52c164`, names the shipped default:

| row | the declared rule (c52c164) | headline before → after | status before → after |
|---|---|---|---|
| `t1_nmc_f2d` | `HEADLINE = "toy"  # the shipped default carries the row's verdict`; config `toy`: "ToyF2() -- the SHIPPED default, i.e. exactly what the generator's `default_inclusive_kernel(deuteron()).tables(x, Q2).f2` returns. THE REPORT ROW'S VERDICT IS THIS CONFIGURATION'S." | toy (8307.9/155) → **ct18nlo** (3401.8/146 on Q² ≥ 1.677) | FAIL → FAIL |
| `t1_nmc_f2d_over_f2p` | the same words ("ToyF2() -- the SHIPPED default; the numerator is exactly `default_inclusive_kernel(deuteron()).tables(x, Q2).f2`. THE REPORT ROW'S VERDICT IS THIS CONFIGURATION'S."), `HEADLINE = "toy"` | toy (2161.4/211) → **ct18nlo** (215.5/189, p = 0.0904) | **FAIL → PASS** |
| `t1_hermes_f2d` | "ToyF2() -- the SHIPPED default … THE REPORT ROW'S VERDICT IS THIS CONFIGURATION'S, ON THE DEUTERON"; `HEADLINE = "toy"  # the shipped default carries the row's verdict ...` | toy (5303.5/68) → **ct18nlo** (1788.2/45) | FAIL → FAIL |
| `t1_g1d_world` | "ToyG1() on ToyF2 -- the SHIPPED default (`--pol-sf toy`, `--unpol-sf toy`). THE REPORT ROW'S VERDICT IS THIS ONE'S." | toy → **nnpdfpol** (45.0/15, 52.6/15, 51.0/28) | FAIL → FAIL |
| `t5_sum_rules` | "the failing sub-row is the shipped toy's Bjorken integral"; the status rule (fail if any evaluated sub-row fails) | the line now LEADS with the NNPDFpol11 Bjorken (0.2084, +12 %, FAIL), the toy's (divergent) recorded after it | FAIL → FAIL |
| `t5_delta_moment_bound` | "(a) the SHIPPED Delta satisfies \|M2\| <= B …, (b) every helicity-state rate of the shipped 6Li kernel is >= -1e-12" | `shipped_kernel()` = `lipolgen.shipped_inclusive_kernel(li6())`; the letter's `moment_A` / `moment_B` keep the toy F₁ through `letter_kernel()` (they reproduce the letter's printed A = −0.2924 on its toy footing, `a_reproduced` still true) | PASS → PASS (\|M₂\|/B 0.0760 / 0.0907 / 0.1023; w_min 0.9135) |
| `t5_weighted_unweighted` | "THE PHYSICS POINT. The shipped defaults of `lipolgen.make_config()` -- … and `default_inclusive_kernel(li6())` (Miller b1, ToyF2, ToyG1), which is what `Pipeline` builds" | the kernel is `shipped_inclusive_kernel(li6())`, what `Pipeline` builds | PASS → PASS (pull +0.359) |
| `t4_pythia_ep_closure` | the bridge sample is "what cli.py does, verbatim": the run's `make_config` and its bridge on `cfg.unpol_sf_obj` | the bridge side on the default (CT18NLO), wired as `cli.py` does | PASS → PASS (opt-in; §10) |
| `t4_djangoh_rad_noRad` (row 7) | the hashed design block (lines 41–187, never edited) needs no LiPolGen Born for its verdict but the recorded sub-rows and the 2026-09-23 history do | the toy PASSED EXPLICITLY at its Born (`unpol_sf="toy", pol_sf="toy"`), a dated D7 addendum after the block | unchanged, bit for bit (`test_bench_djangoh_eltail.py`'s rel 2e−7 pins hold; `test_bench_chain_rc.py`'s "< 0.4" bound not moved) |

**The rule split the check named (its wrong statement 4), recorded, not
resolved:** the T1 rules define the verdict configuration twice — "the
SHIPPED default" and "exactly what `default_inclusive_kernel(deuteron())`
returns". Design A keeps the one-argument kernel the class-level toy, so the
two halves now name different backends. The task's amendment 8 instructs the
"shipped default" reading, and the headlines follow it; the toy rows stay
recorded beside them with their old numbers, so either reading can be read
off the report. The F₂ᵈ/F₂ᵖ flip is robust to the cut: on the frozen shipped backend the plain
Q² ≥ 1 cut gives 244.41/211, p = 0.057 (the verifier, 2026-11-03; the 245.19/211,
p = 0.053 this record quoted was the check's number on the CONTINUED, pre-freeze
backend). The declared verdict drops every point below 1.677 GeV², and
does not cover the frozen region, where 39 % of the default rate now sits.
The author may want the other reading (§13). **The same clause read two
opposite ways** (fix stage, 2026-11-03, physics finding 4): I6 kept Deeps'
verdict on the Hulthén its declaration names (`deuteron_channel()`), reading
"the shipped default" as the declared configuration; I7 re-pointed the T1
headlines to the shipped default after the `ct18nlo` sub-rows were already in
the report. No verdict is wrong as a result (Deeps fails on AV18 too,
579.12/187); which reading the rules mean is the author's.

PLOT_SLOTS and labels follow: the shipped backend takes the "default"
palette slot and the "(shipped since D7)" label in `t1_nmc_f2d`,
`t1_nmc_f2d_over_f2p`, `t1_hermes_f2d`, `t1_g1d_world`, `t5_sum_rules`;
`_plots.py`'s palette comment and usage label, `t1_r_world`'s F₂ labels (its
verdict is R's) and `t2_bonus_spectator_shape`'s species-ratio sub-row
("ToyF2, the shipped species ratio until D7") are labels only.
validation/benchmarks/data/delta_moment_bound_sources.json line 164 ("the shipped
6Li kernel's per-nucleon F1 (ToyF2)") is sha256-pinned and NOT edited; the
harness header says so beside `letter_kernel()`.

## 7. Python: pass / skip / fail per file, before and after (amendment 11)

"The flip alone" is the c52c164 tree with only the two field defaults
flipped and the tier attaching them (83 failed + 17 errors = 100), the state
every re-take in §4 answers. Every file is green after; every cell of the
knob matrix runs.

| file | c52c164 (p / s / f) | the flip alone | after I7 | time c52c164 → after |
|---|---|---|---|---|
| `test_b1_model.py` | 51 / 0 / 0 | 33 / 0 / 18 | 53 / 0 / 0 | 44.0 → 43.4 s |
| `test_bench_boeglin.py` | 13 / 0 / 0 | 13 / 0 / 0 | 13 / 0 / 0 | 1.8 → 1.8 s |
| `test_bench_bonus.py` | 26 / 0 / 0 | 26 / 0 / 0 | 26 / 0 / 0 | 11.6 → 11.5 s |
| `test_bench_chain_rc.py` | 4 / 0 / 0 | 2 / 0 / 2 | 4 / 0 / 0 | 20.4 → 20.4 s |
| `test_bench_deeps.py` | 16 / 0 / 0 | 16 / 0 / 0 | 16 / 0 / 0 | 2.4 → 2.5 s |
| `test_bench_djangoh_eltail.py` | 18 / 1 / 0 | 16 / 1 / 2 | 18 / 1 / 0 | 42.0 → 41.2 s |
| `test_bench_li6_ff_theory.py` | 18 / 0 / 0 | 18 / 0 / 0 | 18 / 0 / 0 | 1.9 → 2.1 s |
| `test_bench_nuclear.py` | 10 / 0 / 0 | 10 / 0 / 0 | 10 / 0 / 0 | 1.7 → 1.7 s |
| `test_bench_report.py` | 25 / 1 / 0 | 25 / 1 / 0 | 25 / 1 / 0 | 3.0 → 3.1 s |
| `test_bench_rfy.py` | 16 / 0 / 0 | 16 / 0 / 0 | 16 / 0 / 0 | 1.0 → 0.9 s |
| `test_bench_spin_deuteron.py` | 24 / 0 / 0 | 24 / 0 / 0 | 24 / 0 / 0 | 13.1 → 12.5 s |
| `test_bench_static.py` | 18 / 0 / 0 | 18 / 0 / 0 | 18 / 0 / 0 | 5.2 → 5.2 s |
| `test_bench_t1_nucleon.py` | 49 / 1 / 0 | 48 / 1 / 1 | 49 / 1 / 0 | 3.6 → 3.6 s |
| `test_bench_t1_r_world.py` | 51 / 0 / 0 | 51 / 0 / 0 | 51 / 0 / 0 | 6.3 → 6.2 s |
| `test_bench_t4_t6.py` | 18 / 1 / 0 | 16 / 1 / 2 | 18 / 1 / 0 | 20.1 → 22.9 s |
| `test_bench_t5_closures.py` | 19 / 0 / 0 | 19 / 0 / 0 | 19 / 0 / 0 | 5.7 → 5.7 s |
| `test_bench_t5_delta_bound.py` | 23 / 0 / 0 | 20 / 0 / 3 | 23 / 0 / 0 | 9.6 → 14.6 s |
| `test_cli.py` | 35 / 0 / 0 | 34 / 0 / 1 | 35 / 0 / 0 | 4.2 → 7.8 s |
| `test_cluster_config.py` | 27 / 0 / 0 | 27 / 0 / 0 | 27 / 0 / 0 | 3.1 → 3.2 s |
| `test_d7_shipped_kernel.py` | — | — | 5 / 0 / 0 | 0.0 → 2.2 s |
| `test_delta_convention.py` | 6 / 0 / 0 | 6 / 0 / 0 | 6 / 0 / 0 | 0.1 → 0.4 s |
| `test_doc_link_gate.py` | 112 / 0 / 0 | 103 / 0 / 9 | 112 / 0 / 0 | 24.4 → 23.7 s |
| `test_hepmc.py` | 6 / 0 / 0 | 6 / 0 / 0 | 6 / 0 / 0 | 0.1 → 0.5 s |
| `test_hfs.py` | 7 / 0 / 0 | 0 / 0 / 7 | 7 / 0 / 0 | 0.6 → 0.8 s |
| `test_kernel_reference.py` | 11 / 0 / 0 | 11 / 0 / 0 | 11 / 0 / 0 | 0.0 → 0.0 s |
| `test_knob_provenance.py` | 517 / 145 / 0 | 495 / 149 / 18 | 517 / 145 / 0 | 200.3 → 361.9 s |
| `test_li6_unpol_band.py` | 12 / 0 / 0 | 12 / 0 / 0 | 12 / 0 / 0 | 18.5 → 17.3 s |
| `test_li7_rank2.py` | 27 / 0 / 0 | 25 / 0 / 2 | 27 / 0 / 0 | 7.6 → 9.0 s |
| `test_meta_provenance.py` | 8 / 0 / 0 | 8 / 0 / 0 | 8 / 0 / 0 | 2.8 → 4.0 s |
| `test_module.py` | 89 / 0 / 0 | 89 / 0 / 0 | 89 / 0 / 0 | 3.4 → 5.1 s |
| `test_mstw_sf.py` | 25 / 0 / 0 | 25 / 0 / 0 | 25 / 0 / 0 | 0.3 → 0.3 s |
| `test_o5_reach.py` | 18 / 0 / 0 | 18 / 0 / 0 | 18 / 0 / 0 | 0.7 → 0.7 s |
| `test_pipeline.py` | 19 / 0 / 0 | 19 / 0 / 0 | 19 / 0 / 0 | 2.5 → 4.7 s |
| `test_pom_set_band.py` | 9 / 0 / 0 | 2 / 0 / 7 | 9 / 0 / 0 | 4.8 → 4.7 s |
| `test_price_defaults.py` | 1 / 0 / 0 | 0 / 0 / 1 | 1 / 0 / 0 | 3.1 → 3.4 s |
| `test_pzz_mode.py` | 25 / 0 / 0 | 22 / 0 / 3 | 25 / 0 / 0 | 1.2 → 4.1 s |
| `test_rc.py` | 32 / 0 / 0 | 32 / 0 / 0 | 32 / 0 / 0 | 58.7 → 64.9 s |
| `test_rc_low_x_anchor.py` | 5 / 0 / 0 | 4 / 0 / 1 | 5 / 0 / 0 | 2.3 → 4.3 s |
| `test_release_metadata.py` | 8 / 0 / 0 | 8 / 0 / 0 | 8 / 0 / 0 | 0.0 → 0.0 s |
| `test_sampler.py` | 9 / 0 / 0 | 9 / 0 / 0 | 9 / 0 / 0 | 0.1 → 0.1 s |
| `test_sf_backend.py` | 54 / 0 / 0 | 31 / 0 / 23 | 57 / 0 / 0 | 9.1 → 10.6 s |
| `test_spdx_headers.py` | 6 / 0 / 0 | 6 / 0 / 0 | 6 / 0 / 0 | 0.0 → 0.0 s |
| `test_throughput.py` | 5 / 0 / 0 | 5 / 0 / 0 | 5 / 0 / 0 | 1.9 → 3.4 s |
| **total** | 1472 / 149 / 0 | 1368 / 153 / 100 | 1482 / 149 / 0 | 544.5 → 737.3 s |

+10 passed: five in `test_d7_shipped_kernel.py`, two D7 tests in
`test_b1_model.py`, three in `test_sf_backend.py`
(`test_the_toy_named_is_the_pre_d7_default_bit_for_bit` on the inclusive,
coherent and tagged channels; the pre-flight, kernel and struck-slot rules
went into the existing G-tests). The 149 skips are the
same tests: the knob matrix's 145 (whose cells were relabelled: 52 ids
renamed — the SF, kernel, `b1_unpol`-on-`li6-convolution`, cdks ×2 and
`inclusive_b1` cells across the specs — every one mapped to its old id with the same outcome, none dropped) and
four opt-in rows.

**The knob matrix, cell by cell** (`i7/knob/m_knob_cells.py`, the matrix's own
`_run` without its assertions, on the c52c164 build with c52c164's test file
and on this tree): 627 cells; base status, variant status, "moved" and the
row's status identical for 626, the relabelled ones mapped by hand
(`ct18nlo -> toy` ↔ `toy -> ct18nlo`, …). The one difference before the last
edit — `tagged-6Li-alpha` `inclusive_b1 off -> on`, whose VARIANT stopped
running at x_max 1.0 (§3.5) and would have passed through the "refused"
branch measuring nothing — now runs at x_max 0.95 on every spec with the same
outcome as on c52c164 (read on the α tag, not-read elsewhere). The cdks ×2
cell runs at 0.85 (it passed through the refused branch at 0.95 under D7, the
check's missed site).

**Wall time.** The full suite 544.5 → 737.3 s (+35 %), essentially all of it
in Pipeline builds on the grids (≈ 0.15 s more each, §3.4): the knob matrix
200 → 362 s. The bench subset: §9.

## 8. C++

| state | cases | assertions | wall |
|---|---|---|---|
| c52c164 | 422 / 422 | 17 244 315 | 480.6 s |
| the flip alone (defaults flipped, no test touched) | 371 / 422 (51 failed) | 40 failed | 455.9 s |
| after I7 | **427 / 427** | **17 244 377** | 444.7 s |

The 51 flip-only failures were every machinery test that builds a
`Pipeline` from a default `PipelineConfig` with no filler installed — the
refusal by name, as designed; `toy_config()` answers them (§4). The +62
assertions are attributed: +59 in the five new D7 cases of
`tests/test_lhapdf.cpp`, +3 in the re-taken D2 case (`-tc` counts on both
builds: T11 1849 = 1849, T16 53 = 53). The build: zero warnings.

## 9. Gates, in the design's order (amendment 11)

| gate | c52c164 | after I7 |
|---|---|---|
| build (`cmake --build build -j2`) | 0 warnings | **0 warnings** (examples included) |
| C++ (`./build/lipolgen_tests`) | 422 / 422, 17 244 315 assertions, 480.6 s | **427 / 427, 17 244 377**, 444.7 s |
| targeted pytest (the re-taken files, run as each was edited) | — | green (§4) |
| full pytest | 1472 / 149 / 0, 544.5 s (wall 545.4 s) | **1482 / 149 / 0, 737.3 s (wall 738.3 s)** |
| bench subset (`LIPOLGEN_BENCH=1 -k bench`) | 351 / 1 / 0, 332.5 s | 350 / 1 / **1** before the report was regenerated (the report test, by design) → **351 / 1 / 0**, 337.7 s after |
| strict citation gate | 1241 / 95 / 7 / 0 / 24; SPIN32 19 / 115 / 0 / 0 / 10; PYTHIA_BRIDGE 0 / 0 / 8 / 0 / 0 | **the same**, after re-reading and re-recording five EDITED blocks (docs/USAGE.md lines 11–66, 163–533 and 1978–2169, include/lipolgen/xsec.hpp lines 197–233, tests/test_pipeline.cpp lines 608–634: each still evidences its citing row); no `--fix` needed |
| R5 (`--audit-ranges`) | 209 on 179 blocks, 0 REFUSED | **209 on 179 blocks, 0 REFUSED** (a note added to `PHYSICS_CHANNELS.md` row 208 made it REFUSED and was withdrawn) |
| R6 (`--records-only`) | 1432 citations in 77 records: 933 / 61 / 110 / 129 / 199, 0 broken | **1432 in 78 records (this one cites in prose only): 924 / 61 / 119 / 129 / 199, 0 broken** |
| SPDX | 133 files | **136** (the three new source files) |
| report `--compare`, `--check` | 13 / 10 / 0 | **14 / 9 / 0**, 34 changes (§10), `--check` OK |

**R6, the record citations D7's edits moved.** Ten broke (content edits on
the same lines, not shifts): nine are annotated `(as of c52c164)` —
run_2026-09-06/phase_D_numbers.md line 1018 (src/core/pipeline.cpp line 1110,
twice) and plan_D7_backend_defaults.md lines 414 (src/lhapdf/lhapdf_sf.cpp lines
39 and 157), 591 (python/README.md 212–213), 659 (t1_nmc_f2d.py 439), 660
(t1_nmc_f2d_over_f2p.py 452), 668 (t5_weighted_unweighted.py 143–145) and
813 (README.md 160–167) — each passing the token bar at c52c164
(`test_every_as_of_anchor_passes_the_bar_at_its_commit`); the census pin in
`python/tests/test_doc_link_gate.py` 110 → **119** with a dated comment naming
D7; `ASOF_READ_BY_HAND` untouched (9). The tenth,
run_2026-09-03/phase_D_numbers.md line 2122, was re-pointed instead (§13).

## 10. The report (`make_report.py --all --plots --compare <pre REPORT.json>`, then `--check`)

**14 pass, 9 fail, 0 blocked** (13 / 10 before) — the one STATUS change is
the expected flip, `t1_nmc_f2d_over_f2p` FAIL → PASS. `--compare`: 34 changes,
every one explained:

* MOVED headlines (the declared rules, §6): `t1_nmc_f2d`, `t1_nmc_f2d_over_f2p`,
  `t1_hermes_f2d`, `t1_g1d_world`, `t5_sum_rules`, `t5_delta_moment_bound`,
  `t5_weighted_unweighted`, and `t4_pythia_ep_closure` (worst 0.076 → **0.065**,
  W1 n_ch_eta[−inf, −2), the plan's §3.8 prediction; ⟨n_ch⟩ 1.009 / 0.994);
* SUBMOVED: the CT18NLO sub-rows' `all` χ² (they include the points below
  Q² = 1.677, now frozen: F₂ᵈ 3523.99 → 3959.37, F₂ᵈ/F₂ᵖ 295.36 → 293.00,
  g₁ᵈ COMPASS / HERMES / E143 / SMC 52.6 → 62.9, 78.4 → 79.7, 64.61 → 64.56,
  15.2 → 27.9 on the ToyG1-on-CT18NLO row, NNPDFpol SMC 13.42 → 13.69 — the
  DIS-cut numbers do not move), the T4 sub-rows, the t5 positivity sub-row
  (w_min 0.9494 → 0.9135, A_T [−0.173, 0.0865]);
* `t6_hepmc3_convention` C04 `n_checked` 1370 → 1369 (the number of nuclear
  PDG codes in its 8 `lipolgen-run` files — the default event stream moved;
  0 violations, PASS).

Unchanged: rows 7 (`t4_djangoh_rad_noRad`, its figure byte-for-byte but the
footer), 10, 11, 12, 16 and every T2 / T3 row. `--check` **OK**. Figures: the
seven whose content moved are kept (`t1_g1d_world`, `t1_hermes_f2d`,
`t1_nmc_f2d`, `t1_nmc_f2d_over_f2p`, `t4_pythia_ep_closure`,
`t5_delta_moment_bound`, `t5_sum_rules`); twelve whose only change was the
footer stamp were restored. `price_defaults.png` regenerated once from the
full table (N = 200 000, seed 20261026, six families, 358.6 s — the
`li6-incl-r` family 325.9 s of it: `li6-convolution` on `--b1-unpol ct18nlo`
evaluates CT18NLO inside the convolution); the li6-incl-r band ×2 row runs at
x_max 0.85. README's report sentence reads 14 / 9.

## 11. Documents and CI (amendments 9, 10)

Every edit to a cited file is **line-neutral** (old and new text with the
same number of lines, checked per edit) except where nothing below the edit
is cited; new text went to the ends of sections with nothing cited after
them.

* **Code comments and docstrings**: `pipeline.hpp` (the `UnpolSfSource` /
  `PolSfSource` / `B1UnpolSource` blocks, the field and `kernel` comments, the
  `StruckClusterOptions` and `default_inclusive_kernel` docs — "default" →
  "class-level" where the null branch is meant, the grid clause on the
  freeze, the "neutron-tagged" run corrected to the struck-neutron
  `tagged-d-p`, the 7Li 590952.42641509 labelled toy beside 494062.76099942),
  `pipeline.cpp`, `sf.hpp`, `xsec.hpp`, `pythia_bridge.hpp`,
  `b1_nuclear.hpp`, `rc.hpp`, `lhapdf_sf.cpp` (header: the shipped defaults,
  LO, the freeze), the bindings' enum and field docstrings,
  `lipolgen/__init__.py` (`UNPOL_SF` / `POL_SF` / `B1_UNPOL`, `make_config`,
  `run`, the `make_plan` ValueError's measured zero), `cli.py` (every help that
  called the toy the default; `--x-max` names the six refusals, the ×2 row and
  `--inclusive-b1`; the `--inclusive-b1` A_zz prices re-taken; the tier
  refusals; `B1_UNPOL_GATE_ROW`; `b1_gate_lines`), `CMakeLists.txt` (the
  option text and the missing-`lhapdf-config` WARNING: every default run of
  such a build is refused by name).
* **`docs/USAGE.md`**: the C++ examples at 140–160, 224–235, 1279–1290 and
  1713–1725 include `lhapdf_sf.hpp` and install / attach the backends (the
  `li6-convolution` one names `--b1-unpol ct18nlo`); compiled and run as
  written (`i7/usage_cxx/ex.cpp`: 496334.155954, 496334.155133 with the b₁
  F₁ the same object, 9424.476321, the toy 591846.168878). 282–300 (the
  denominator's F₁), §2a's gate and b₁-unpol table, the top-x-cell figures,
  §2b re-titled and its price table re-taken on D7 (tagged rows on the VMC
  wave), the grid clause on the freeze (39.23 %, 48.37 %, 496334.2 → 301599.8
  pb at q2_min 1.7), "What is refused" (the lazy grid slots, the kernel and
  struck-slot rules, `li6-convolution`), the ⁷Li 590952.42641509 labelled toy,
  the knob-provenance banner re-run verbatim (`--pol-sf toy`; 69 knobs —
  the page's 68 was stale before D7), the I5 / I6 numbers of §5, and a NEW
  **§10** (appended): the defaults, LO, the freeze and its `meta`, the
  HepMC3 gap, prices, the refusal table, the deuteron split, data and the toy
  escape, the C++ API contract, what pins it.
* `README.md` (39–41 the reference JSONs pin the class-level kernels, the
  selector table, 160–167, 186–198, 218–222, 243–256); `docs/PHYSICS_CHANNELS.md`
  (12–13, 120, 128, 136–143, 191, 196, 198–199, the I5 RC rows 384–385 and
  828–831); `docs/PACKAGING.md` (the set pins in §1.1, the wheel table's three
  rows, a NEW §2.5: the data half gates the default, no `lhapdf.conf` moves
  it, the pins, and **the C++ API contract**); `docs/OPEN_ITEMS_SOLUTIONS.md`
  (board row 14 and the table row 92, §10's "shipped default ToyF2" → `--b1-unpol`'s
  default, §9's RC numbers with D7 beside the toy, §14's heading, the
  tagged A_zz and the ⁷Li zero); `docs/benchmarking/08_channel_priority.md`
  (36–39, 69, 206, 233, 253, 327, 333, 361–364 and a D7 log row);
  `BENCHMARK_PLAN.md` (a dated D7 status sentence); 00_in_tree_checks.md line 244,
  02_data_nucleon_deuteron.md line 277, docs/references/01_nucleon-deuteron-data.md line 472;
  `validation/benchmarks/README.md` (F-1, F-2 PASS, D-2/D-1/D-6, F-7, T5 rows);
  `python/README.md`; `docs/T2_CHAIN.md` (the T1 numbers on D7).
  `docs/HEPMC3_CONVENTION.md`, `validation/reference/*.json`,
  `docs/theory/SPIN32_FINITE_GAMMA.md` and the hashed design block of
  `t4_djangoh_rad_noRad.py`: not edited.
* **The registry** (`docs/open_items/run_2026-09-03/STATUS.md` rows 4 and
  20, and the board rows 14 / 15): the dated registry is history and is not
  edited here; its D7 texts are proposed in §14.
* **CI**: `.github/scripts/build_deps.sh` verifies the three set tarballs'
  SHA-256 (`declare -A SET_SHA256`; re-enacted: the three hash ok, a
  truncated `CT18NLO.tar.gz` fails with "checksum mismatch", exit 1);
  `.github/workflows/ci.yml` gains the step "the shipped LHAPDF sets absent —
  refused by name, never the toy (D7)" in `test` (re-enacted line for line on
  this machine: exit 1 with LHAPDF's text and `--unpol-sf toy` in the
  message under both an empty data path and a conf-only one, the toy run
  exits 0, pytest returns 3), and the wheel smoke runs the toy NAMED with
  nothing set (`<x>` = 0.027023337580082865, the pre-D7 number to the bit),
  the default refused under an empty data path (caught, exit 0) and the
  default on the restored sets (`<x>` = 0.033169118940819454, `meta` ct18nlo)
  — re-enacted with the in-tree build under `env -i`; the auditwheel wheel
  itself was not rebuilt.

## 12. The check's items: applied, and not

**Its twelve amendments, as the task restated them — all applied:** (1) the
per-PDF setters, the clamp, the opposite-`lhapdf.conf` test and the
extrapolation / below-grid fractions (both backends) in `meta`; (2) the
pre-flight keeps "Info file not found for PDF set", runs before any LHAPDF
`Config`, names the set and the escape; (3) one kernel-beside-selectors
rule, struck slots `Custom`-only, the three reason strings rewritten; (4) the
lazy attach, nothing loaded by `PipelineConfig()` / `make_config()` / the CLI
import, the timing table (§3.4); (5) `test_d7_shipped_kernel.py`; (6) the
refusals documented (USAGE §10, the `--x-max` help, §3.5) and the positivity
text computing its cure, `price_defaults`' `li6-incl-r` naming
`--b1-unpol ct18nlo` and its ×2 band row at 0.85; (7) row 7 on the toy at its
Born, its pins and the chain_rc bound untouched; (8) the verdict rules
followed and quoted (§6), `letter_kernel()` keeps the toy F₁, PLOT_SLOTS and
labels, the I4 / I5 / I6 numbers re-measured (§5); (9) the text sites (§11)
and the CI pins, the sets-absent step and the wheel step, re-enacted; (10)
provenance in `meta` and the banner, the HepMC3 gap stated, the C++ API
contract in USAGE §10 and PACKAGING §2.5; (11) the gates in order with the
per-file tables and the wall time (§7–§9); (12) this record.

**Its 21 missed sites — applied**, each where §11 lists it, with three
qualifications: delta_moment_bound_sources.json line 164 is left (sha256-pinned)
and noted in the harness header; registry row 20 (and row 4) of the dated
2026-09-03 registry are not edited — §14 proposes their text;
`validation/b1_li6_table.py`'s 216 / 277–286 / 380 keep "the ToyF2 default
row": that row IS `--b1-unpol`'s default, which D7 kept (its header's
"SHIPPED DEFAULT" sentence at 30–34 was reworded).

**Its 9 wrong statements — each answered:** (1) the reason strings; (2) the
rule is one rule, and the caller-kernel tests sit beside the SHIPPED
selectors (T11 / T16, `test_b1_model.py`'s kernel sites) rather than the
toy; (3) the toy at row 7's Born, so neither the history nor the 0.4 bound
moves; (4) recorded as an open point (§6, §13); (5) the flip-state counts are
in §7; (6) the deuteron split measured and the development-item text
corrected (§1, §3.5); (7) `--inclusive-b1` on the ⁶Li α tag documented
everywhere the other refusals are; (8) the default no longer depends on the
environment (§2.3); (9) the toy `make_config` re-measured: 0.06–0.07 ms.

**Its 10 underrated risks — each closed or stated:** 1, 8 (§2.3); 2 (§2.4,
`lhapdf_quiet` guarded); 3 (`test_d7_shipped_kernel.py`, and USAGE §10 says
what the reference JSONs pin); 4 (the cdks ×2 cell runs at 0.85; G1b and
`cli.DEFAULTS` asserted by literals; the wheel step reworked); 5 (row 7 on
the toy Born); 6 (nine record citations annotated `(as of c52c164)`, one
re-pointed, the census re-pinned — §9); 7 (lazy attach: no best-effort load,
no enum-assignment trap); 8 (§3.6); 9 (the sampler's message); 10 (stated as
a consequence of D7 in USAGE §10, PACKAGING §2.5 and `conftest.py`).

**Not applied:** the plan's OPTIONAL `meta` record of the loaded member's
identity (`unpol_sf_set`: name, member, `DataVersion`, member-file sha256) —
the tarballs are pinned in `build_deps.sh` instead, and the run records the
set name and the extrapolation; the plan's "stricter alternative" of a
configure-time `FATAL_ERROR` when `lhapdf-config` is missing — the WARNING
is reworded instead (a build-surface change for the author, §13); the
design's `--fix` pass was not needed (no citation MOVED in the strict
documents: every code and document edit was line-neutral).

## 13. Open points, for the author and the supervisor

1. **The T1 rules' split** (§6): the headlines follow "the shipped default";
   the toy rows stay recorded. If the author reads the rules as "the
   one-argument `default_inclusive_kernel`", F₂ᵈ/F₂ᵖ is FAIL again and the
   report 13 / 10.
2. **`t5_delta_moment_bound`'s `m2_numerical_error < 1e-15`**: measured on
   the analytic toy F₁ (~1e−19), it is kept on that integrand
   (`letter_kernel()`, the harness's own rule); the shipped row's error on
   CT18NLO's C¹ interpolant (3.3e−14 / 4.9e−14 / 5.5e−14, 2–4e−10 relative)
   was (I7) held to the harness's DECLARED `QUAD_REL_TOL` = 1e−7. That IS a test
   bound moved after a number was seen (I7 wrote "No bound was moved after a
   number was seen", which was false for this assertion; corrected by F7,
   §16, V7-5): the shipped row's pin was `e < 1e-15` absolute before D7, the
   CT18NLO values fail it, and the declared convergence rule it then carried
   is about 10⁴ looser in absolute terms. F7 (2026-11-03) removed it: the grid
   row's error is now RECORDED with no pass bound; any bound is the author's.
3. **A missing `lhapdf-config`**: the CMake WARNING now says every default
   run of such a build is refused; the plan's stricter alternative, a
   `FATAL_ERROR` unless `-DLIPOLGEN_WITH_LHAPDF=OFF` is explicit, is the
   author's call (it ends the "optional tier" convention at configure time).
4. **The record citation re-pointed, not annotated**:
   run_2026-09-03/phase_D_numbers.md line 2122 cited pipeline.hpp lines 628–635, the
   grid clause, for the "flag pair reaching EVERY kernel" sentence, which is
   at lines 694–696 before and after D7 (it passed R6 on the single token
   "flag"); it now points there. Nine other record citations D7's edits moved
   are annotated `(as of c52c164)` (§9).
5. **`b1_unpol` is eager, not lazy** (§2.2): a deliberate departure from
   "every grid slot attaches lazily", to keep the D2 contract that an empty
   `Ct18Nlo` b₁ slot is refused.
6. **HepMC3 carries no record** of the structure functions (nor of D4's
   convention, D5's tail model, D6's wave): one run-level attribute would
   change `docs/HEPMC3_CONVENTION.md`, which T6 pins by sha256.
7. **EPPS21 and R1998** stay out of D7, pending the author.
8. **CI**: the new step and the reworked wheel step are re-enacted locally,
   never run on GitHub (registry row 14); the auditwheel wheel was not rebuilt.

## 14. Registry and board texts — proposed, not applied

* `run_2026-09-03/STATUS.md` **row 4** (the ban lift), append: *(D7,
  2026-10-27: CT18NLO became the shipped `--unpol-sf`. The 0.440 row is
  `--b1-unpol`'s default `toy`, which D7 kept; under the shipped `ct18nlo`
  the `li6-convolution` branch must name `--b1-unpol` — `ct18nlo` (G3b
  0.719, inside the window but not CDKS's PDF) or `mstw` (0.843, the lift's
  configuration) — so no default run sits on the 0.440 row any more.)*
* **row 20** (§15.5 D1–D13), append: *(D7: "the shipped `toy` default" of
  §15.5 D2 is `--b1-unpol`'s, which D7 kept; the RATE's default is CT18NLO
  since 2026-10-27.)*
* **row 15** (publish the wheel), append: *(D7 consequence: a pip-installed
  default run needs the CT18NLO / NNPDFpol11_100 data — refused by name
  without them; the toy named runs with no data. Deciding to publish stays
  the author's.)*
* `DECISIONS.md` D7: the decision text stands as recorded; the
  development-item wording of §1 (with the measured deuteron split) is
  offered as its gloss.

## 15. The files

**New** (4 + this record): `src/lhapdf/shipped_sf.cpp` (the set names, the
two cached builders, `attach_shipped_sf_backends`,
`install_lhapdf_sf_filler`), `src/pythia/bridge_config.cpp`
(`check_bridge_species_backend`, `set_pythia_hadronizer`'s C++ rule),
`tests/sf_test_support.hpp` (`toy_config()` for the core suite),
`python/tests/test_d7_shipped_kernel.py`.

**Code**: `include/lipolgen/{pipeline.hpp, lhapdf_sf.hpp, pythia_bridge.hpp,
sf.hpp, xsec.hpp, b1_nuclear.hpp, rc.hpp}` (the first three carry the D7 API,
appended after every cited line; the rest comments only), `src/core/pipeline.cpp`
(`validate()` line-neutral, the ctor's fill, the registry and `toy_config`
appended), `src/core/sampler.cpp` (the refusal text, line-neutral),
`src/lhapdf/lhapdf_sf.cpp` (rewritten around the cited declarations, which
keep their lines), `python/bindings.cpp`, `python/lipolgen/{__init__.py,
cli.py}`, `CMakeLists.txt`, the four `examples/generate_*.cpp`.

**Tests**: `tests/{test_lhapdf, test_pipeline, test_b1_nuclear,
test_rc_pipeline, test_t2, test_breakup, test_fsi, test_triton_sf}.cpp`;
`python/tests/{conftest, test_sf_backend, test_b1_model, test_knob_provenance,
test_li7_rank2, test_pzz_mode, test_rc_low_x_anchor, test_price_defaults,
test_cli, test_hfs, test_pom_set_band, test_bench_t1_nucleon,
test_bench_t5_delta_bound, test_bench_t5_closures, test_bench_bonus,
test_doc_link_gate, test_li6_unpol_band}.py`.

**Validation**: `validation/price_defaults.py`, `validation/b1_li6_table.py`,
`validation/physics_channels_ranges.json` (five re-recorded fingerprints),
`validation/benchmarks/{t1_nmc_f2d, t1_nmc_f2d_over_f2p, t1_hermes_f2d,
t1_g1d_world, t1_r_world, t2_bonus_spectator_shape, t4_djangoh_rad_noRad
(after its hashed block only), t4_pythia_ep_closure, t5_delta_moment_bound,
t5_sum_rules, t5_weighted_unweighted, _plots}.py`, `validation/benchmarks/README.md`.

**Documents**: `README.md`, `python/README.md`, `docs/{USAGE, PHYSICS_CHANNELS,
PACKAGING, OPEN_ITEMS_SOLUTIONS, T2_CHAIN}.md`,
`docs/benchmarking/{00_in_tree_checks, 02_data_nucleon_deuteron,
08_channel_priority, BENCHMARK_PLAN}.md`, `docs/references/01_nucleon-deuteron-data.md`;
the regenerated `docs/benchmarking/REPORT.{md,json}` and eight figures
(`price_defaults`, `t1_g1d_world`, `t1_hermes_f2d`, `t1_nmc_f2d`,
`t1_nmc_f2d_over_f2p`, `t4_pythia_ep_closure`, `t5_delta_moment_bound`,
`t5_sum_rules`); three dated records' citations (§9); `STATUS.md`.

**CI**: `.github/scripts/build_deps.sh`, `.github/workflows/ci.yml`.

**Not touched** (the ground rules): `docs/HEPMC3_CONVENTION.md`,
`validation/reference/*.json`, `docs/theory/SPIN32_FINITE_GAMMA.md`, the
hashed design block of `t4_djangoh_rad_noRad.py`, every vendored data file,
`ASOF_READ_BY_HAND`, the sibling repository.

## 16. Verification (2026-10-30 / 2026-11-03)

The adversarial verification of the tree I7 left (session scratch `i7/verify_result.json`; its evidence directory `v7/` is gone). Its verdict:
I7 "holds up on physics, toy-path bit identity, the freeze, the pre-flight, lazy attach, the knob matrix, every gate tally and the
report". It found one medium finding and nine low or info findings. No default physics number, tolerance, verdict rule, reference
JSON or hashed block is wrong. The fix stage F7 started on 2026-10-30 at about 04:40 UTC. A usage limit stopped it at about
05:40 UTC, while it was running its final gates. It had applied most findings in place, but its scratch directory did not survive
and it wrote no record. The edits it left carry `F7` marks in 20 files, and other F7 edits are recognisable by date (files
modified 2026-10-30 05:00–05:19 UTC). F7 resumed on 2026-11-03. That work established the state of each finding against the tree,
re-checked each one, did what was left, and ran every gate on the final tree. No finding was refuted. Committed as `7fd441c`
(this sentence said "Nothing is committed" until the fix stage of 2026-11-03).

| # | lens | finding | state found in the tree (2026-11-03) | re-check | what was done |
|---|---|---|---|---|---|
| V7-1 (medium) | environment contract | With CT18NLO or NNPDFpol11_100 missing, the four T1 harnesses crashed with KeyError `'dis'` / `'experiments'` (the report row became `error`). `price_defaults.py` crashed on `'fingerprint_tt'`. This contradicts §2.4's "keep working". | **Already applied** by the interrupted F7: each harness's `_headline()` raises LHAPDF's own "Info file not found for PDF set '<set>'" when the verdict configuration is blocked, and a tier-less build is refused by name. `price_defaults` marks every probe of an `environment` default row `environment`. Tests: `test_bench_t1_nucleon.py`'s `test_lhapdf_configs_are_reported_not_dropped` (re-taken), `test_a_missing_verdict_set_blocks_the_row_by_the_environment` (subprocess, `LHAPDF_DATA_PATH=<empty>::`), `test_the_verdict_configuration_without_the_tier_is_refused_by_name`, `test_bench_t1_r_world.py::test_hermes_ct18nlo_reported_not_dropped`, and `test_price_defaults.py`'s hidden-set test. §2.4 is corrected in place. | Confirmed on 2026-11-03 with `make_report.py --only t1_nmc_f2d --only t1_g1d_world` and an empty data path: exit 0, and both rows are `blocked (environment)` ("environment: LHAPDF set CT18NLO / NNPDFpol11_100 not installed"). | **Applied now:** an end-to-end test, `test_make_report_records_the_four_t1_rows_blocked_without_the_sets` (the end of `test_bench_t1_nucleon.py`). It runs make_report on all four T1 rows with the sets hidden, then asserts exit 0 and `blocked (environment)` for each row, naming the right set. **Deviation from the brief's wording:** the HARNESS itself exits **2** in this case, not 0. That is the project's own rule for a row whose verdict configuration needs a missing set (validation/benchmarks/README.md convention 6; rows 4 and 13 since run 2026-10-26). make_report, which the brief's "blocked (environment) row" refers to, exits 0. A stale pin was also found and fixed: `test_hermes_verdict_follows_from_the_numbers` still compared the row status with the `toy` configuration. Both are FAIL, so it passed by coincidence. It now reads `configs_detail[hmod.HEADLINE]` (line-neutral). |
| V7-2 (low) | docs | python/README.md's T2 example is refused on a default config, because the bridge's own ToyF2 does not match `unpol_sf = ct18nlo`. | **Already applied**: `fill_grid_sf_backends(cfg_pipeline)` and `o.f2_source = cfg_pipeline.unpol_sf_obj` (line-neutral). | Read; this is the verifier's own recipe. | Nothing further. |
| V7-3 (low) | provenance | meta and the banner credited a polarised grid (`nearest`, Q² = 1, 0 %) on runs where g₁ was not read. | **Already applied**: `sf_grid_extras(p, knob("pol_sf"))` writes the unread label and NaN where g₁ did not run, and the banner's polarised line needs `pol_read`. `test_sf_backend.py` is re-taken under F7 comments. | The full suite passes, including test_sf_backend.py. | Nothing further (§2.3 updated by the interrupted F7). |
| V7-4 (low) | honesty | The LO and free-nucleon caveats were missing where "shipped" is said: the four T1 headlines, README.md:57, the `--unpol-sf` help (cli.py:622) and pipeline.hpp:587. | **Already applied**: the harness labels read "(shipped since D7, LO)" and the report was regenerated (2026-10-30 05:18 UTC). README, the CLI help, pipeline.hpp and the bindings name "the FREE nucleon (no EPPS21)". | `make_report --check` is OK on the final tree. The four T1 generator strings differ from I7's report only by the added ", LO"; every χ² / p is as the verifier reproduced (3401.8/146, 215.5/189, 1788.2/45, 45.0 / 52.6 / 51.0). | Nothing further. |
| V7-5 (low) | test bounds | `test_bench_t5_delta_bound.py:361`: the shipped row's pin `e < 1e-15` was replaced by `QUAD_REL_TOL·M2` after CT18NLO's 3.3e−14 / 4.9e−14 / 5.5e−14 failed it, and the record said no bound had moved. | **Partly applied**: the record's sentences were reworded (§1, §4, §13 item 2), but the test still carried the bound that had been chosen after the value. | Confirmed against `c52c164`: the old line was `all(e < 1e-15 ...)`. The `QUAD_REL_TOL` rule dates from before D7, but it is the harness's run-time rule for a broken harness, not this test's bound. | **Applied now**, line-neutral: the declared 1e−15 stays, unchanged, on the analytic toy integrand (`_toy_m2_errors_below(mod, 1e-15)` on `letter_kernel()`). The grid row's error is **recorded, without a pass bound** (asserted only to be three finite, non-negative values). The reason is written in the comment above the helper: any bound chosen now would be chosen after the value. §1, §4 and §13 item 2 now say so. |
| V7-6 (low) | docs | USAGE §8's throughput table (2.83 / 1.02 / 2.32 M ev/s) was toy and pre-D7, and was not labelled. | **Already applied**: the table's header reads "toy SFs, 2026-08-30". The line above gives a re-measurement dated 2026-10-30 on D7's grids: inclusive 0.65, tagged ⁶Li α 0.36, coherent 1.08 M ev/s, against 1.92 / 0.585 / 1.13 on `c52c164` (×0.34 / ×0.61 / ×0.96). | The raw timings did not survive (the scratch is gone). The ratios agree with the verifier's independent 0.494 vs 1.62 M ev/s (×0.30, loaded machine). | Nothing further. The 2026-10-30 timings are not re-measured here: the machine was not idle while the gates ran. |
| V7-7 (low) | docs | Stale wording: USAGE.md 950-952 and 1080, USAGE §10's ⟨Q²⟩, t5_delta_moment_bound.py:66, and T2_CHAIN.md:113's unlabelled 19 / 11 / 1. | **Already applied**: "toy g₁ⁿ (shipped until D7)" and "struck-neutron"; ⟨Q²⟩ 4.18 → 4.54 (inclusive) / 4.47 (tagged); t5 line 66 names `shipped_inclusive_kernel(li6())`; T2_CHAIN gives 19 / 11 / 1 on toy SFs and 20 / 11 / 1 on D7's grids. plan_D7's two quotations of t5:66-73 are annotated `(as of c52c164)`, and the census was re-pinned 119 → 121 with a dated F7 comment (test_doc_link_gate.py). | Read at each site. | Nothing further. |
| V7-8 (low) | one place | The four harnesses' `HEADLINE` literals were not tied to `UNPOL_SF_DEFAULT` / `POL_SF_DEFAULT`. | **Already applied**: `_headline()` checks `HEADLINE == unpol_sf_name(UNPOL_SF_DEFAULT)` (`pol_sf_name(POL_SF_DEFAULT)` for g₁ᵈ) and stops the row by name otherwise. `test_headline_is_the_libraries_shipped_default_or_the_row_stops` asserts both the equality and the refusal. | Read; the test passes. | Nothing further. The literal is kept, with an assertion rather than a derivation, because its labels and docstrings name the set. |
| V7-9 (info) | supply chain | `build_deps.sh`'s SET_SHA256 pins verify only on a fetch. | **Already applied**: a set already in the prefix is re-verified through its cached tarball in `$SRC`. Without that tarball the script says "NOT re-verified". The comment states that a CI cache hit skips the script. | Read. | Recorded. Nothing further. |
| V7-10 (info) | message | With the sets present and only `lhapdf.conf` missing, the pre-flight still advised "install the set". | **Already applied**: the cure now follows what is missing ("Put LHAPDF's own lhapdf.conf … in a directory on LHAPDF_DATA_PATH -- the set itself is there"), pinned in `test_d7_shipped_kernel.py`. | The suite passes. | Recorded. Nothing further. |

**Gates on the final tree** (2026-11-03, one at a time):

- **cmake build**: nothing to rebuild, since the build tree post-dates the last C++ edit (the 2026-10-30 build), 0.4 s.
- **`./build/lipolgen_tests`**: 427 cases, 427 passed, 0 failed, 0 skipped, **17 244 377 assertions**, 462.4 s wall. This is I7's tally; F7's C++ edits are messages and meta keys, with no new case.
- **The full pytest**: **1497 passed, 149 skipped, 0 failed, 853.8 s** (wall 854.7 s). That is I7's 1482 plus the 15 tests F7 added: four T1 missing-set subprocess tests, four tier-refusal tests, four `HEADLINE` tests, the make_report end-to-end test, one `price_defaults` hidden-set test and one `lhapdf.conf`-only test.
- **The bench subset** (`LIPOLGEN_BENCH=1 … -k bench`): **364 passed, 1 skipped, 0 failed, 340.6 s**. That is I7's 351 plus the 13 new tests whose ids carry `bench`.
- **Citation gates**: strict **1241 / 95 / 7 / 0 / 24** (SPIN32 19 / 115 / 0 / 0 / 10; PYTHIA_BRIDGE 0 / 0 / 8 / 0 / 0). `--audit-ranges` **209 on 179 blocks, 0 REFUSED**. `--records-only` **1432 citations in 78 records, 922 / 61 / 121 / 129 / 199, 0 broken**: two citations moved from live to annotated, plan_D7's t5:66-73 (census 121). Neither `--fix` nor `--record-ranges` was needed, because every F7 edit in a cited file is line-neutral and the additions sit at file ends.
- **SPDX**: 136 files.
- **`make_report.py --check`**: OK (77.2 s).

**The report.** 23 harnesses: **14 pass, 9 fail**, 0 blocked, 0 blocked (environment), 0 error. This is I7's summary, with the same verdicts and numbers. Against I7's report the only changes are the four T1 generator strings, which gained "LO" (V7-4). The interrupted F7 regenerated the report on 2026-10-30; `--check` confirms it on the final tree. The 2026-11-03 edits touch no harness and no default, so no further regeneration was run. I7's pre-F7 REPORT.json did not survive, so `--compare` against it could not be repeated. With the sets hidden, make_report records the four T1 rows `blocked (environment)` (V7-1).

**Corrected in place** (each with a dated F7 note): §1's sentence on the moved tolerance, §2.3's polarised meta keys, §2.4's "keep working", §4's t5 row, and §13 item 2.

**Left for the supervisor and the author:**

1. Whether the grid row of `t5_delta_moment_bound` should get a declared numerical-error bound at all (none is set now).
2. The brief's "exit 0" for a harness with a missing set. F7 follows the README's convention 6 instead (the harness exits 2; make_report exits 0 with `blocked (environment)`).
3. USAGE §8's D7 throughput, re-measured 2026-10-30, has no surviving raw log; it can be re-taken on an idle machine.
4. §13's open points 1 and 3–8 are unchanged.

**After the run's verification (2026-11-03).** The four-lens verification of
the whole run and its fix stage are recorded in `phase_V_verification.md`.
Two of its items concern this phase directly: (a) a non-null `unpol_sf_obj` /
`pol_sf_obj` under a grid selector is now refused unless it is the tier's
shipped object (code finding 1); (b) the regression lens built
`-DLIPOLGEN_WITH_LHAPDF=OFF` after D7 (Python off, in its scratch): 0
warnings, the examples not built (as `docs/PACKAGING.md` says),
`lipolgen_tests` 411 / 411 with 17 242 996 assertions — no record had
measured it.
