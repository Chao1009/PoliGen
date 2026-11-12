<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-11-10 — phase P2b: the deuteron's b₁ slot confined to the spin-1 positivity domain (J1 and the supervisor's note of 2026-11-12)

Worker P2b, 2026-11-12, on `claude/run-2026-10-26` at `d86fbf2` (phase P2
committed; this phase uncommitted, the supervising session commits). What it
implements: the **"Supervisor's note, 2026-11-12" under J1** in `DECISIONS.md`
— **the SUPERVISOR's choice** (J1 itself is the author's decision; the
author's words on item 1, "go with your recommendation", say nothing about a
treatment). Every number below was measured in this phase on this machine
(shared with two other sessions since 2026-11-12) with `source env.sh`;
scripts and raw logs are in the session scratch `run_1110/p2b/`
(`scripts/m1_dump.py`, `n3.py`, `n45.py`, `n7.py`, `m6_hermes.py`,
`ab_runs.sh`, `ab_compare.py`; the scratch is volatile, every number a
reader needs is here). **pre** = `d86fbf2`'s build, every pre measurement
taken before the first edit; **post** = this tree. Dated comments of this
phase: `J1 note (2026-11-12)`.

## 0. What changed, in one paragraph

ONE kernel option, `InclusiveKernel::Options::confine_spin1_rank2` (include/
lipolgen/xsec.hpp; DEFAULT OFF, so a kernel a caller builds is unconfined and
P2's reading stays constructible from C++ and Python), confines the kernel's
FINAL b₁ per nucleon in `InclusiveKernel::tables` to the leading-order spin-1
positivity domain, **b₁ ← min(max(b₁, −(3/2)F₁), 3(F₁ − |g₁|))**, with the
kernel's own F₁ and g₁ at the same (x, Q²) — g₁ as the vector term uses it
for the fully polarised ion, i.e. the nucleus' g₁ per nucleon after A3's
clip (|g₁|/F₁ = 0.9136 = the deuteron's eff_pol where both nucleons are
clipped, the refusal message's "A1 = g1/F1 is 0.9136"); b₂ in step (an unset
b₂ — the default — stays 2x·b₁, recomputed from the confined b₁; a set b₂ is
scaled by the same factor); an empty window (3(F₁ − |g₁|) < −(3/2)F₁) or
F₁ ≤ 0 leaves b₁ alone for the sampler's refusal to speak; b₃, b₄ (empty →
0 on the default kernel) and Δ (`toy_delta_gluon`, 1e−2) are not touched; a
non-spin-1 ion with the option set throws. `default_inclusive_kernel` sets it
in the DEUTERON branch only (the class toy and the shipped kernel alike).
Visibility mirrors A3: `b1_rank2_confinement(Pipeline)` (pipeline.hpp
`B1RankTwoConfinement`, the definition written there before the measurement)
→ `meta["b1_clipped_rate_fraction"]`, written only where the confinement is
applied (absent elsewhere), one banner line, `meta["rank2_input"]` a clause;
g₁'s share is not in `run_provenance` / HepMC3, so neither is this one. One
consequence the brief did not foresee and the tree's own knob gate forced: g₁
now enters the deuteron's rank-2 slot, so `--pol-sf` is READ on every
inclusive deuteron run, the tensor plans included (`pol_sf_is_read`; §1 N2,
§5 item 2). Measured: all 24 default deuteron plan × config runs build at
x_max 1.0 (0 of 24 on P2's tree); the confinement binds on 40 / 59 / 59 cells,
4.792e−7 / 3.871e−7 / 3.536e−7 of the rate, all on its upper (joint) half;
⁶Li, ⁷Li, tagged, coherent byte-identical; HERMES unchanged to every digit.
The sampler's refusal is unchanged (the backstop).

## 1. Every number or label that moved, before → after

### The supervisor's measurement (2026-11-12, on P2's tree; copied from the supervisor's `joint.json`, labelled as the supervisor's)

Cells of the full x_max 1.0 window outside the single-function bound (S:
b₁/F₁ ∉ [−3/2, 3]) and outside the joint domain (J), and the plan states each
candidate clip would leave refused (CLI fill):

| ion, cfg | cells | max b₁/F₁ | max \|g₁\|/F₁ + b₁/(3F₁) | S: cells, share | J: cells, share, rows | refused state-cells: none / S / J (helicity-flip; tensor plans) |
|---|---|---|---|---|---|---|
| d, 0 | 1950 | 5.4221 | 2.7210 | 15, 1.677e−9 | 40, 4.792e−7, 0.871 + 0.955 | 60 / 32 / 0; 90·15·30 / 0 / 0 |
| d, 1 | 3053 | 6.9931 | 3.2446 | 24, 1.658e−9 | 59, 3.871e−7, 0.871 + 0.955 | 96 / 44 / 0; 144·24·48 / 0 / 0 |
| d, 2 | 3498 | 6.9931 | 3.2446 | 24, 1.563e−9 | 59, 3.536e−7, 0.871 + 0.955 | 96 / 34 / 0; 144·24·48 / 0 / 0 |
| ⁶Li, 0 | 1947 | 1.7759 | 0.8746 | 0 | 0 | 0 / 0 / 0 |
| ⁶Li, 1 / 2 | 3051 / 3498 | 2.2904 | 1.0461 | 0 | 3, 2.372e−11 / 2.712e−11, 0.955 | 0 / 0 / 0 |

### N1 — kernel tables (`m1_dump.py`; P2's M1 points, 8 x × 4 Q²)

Deuteron, shipped kernel (the `Pipeline`'s own kernel gives the same cells
bit for bit): b₁/F₁ pre (|g₁|/F₁; upper bound 3(1 − |g₁|/F₁)) → post where it
binds; "up" = the upper (joint) half; the lower half binds nowhere.

| x | Q² = 2.5 | 5 | 100 | 1897 |
|---|---|---|---|---|
| 0.012 – 0.7 | inside (e.g. 0.7: 0.03756 (0.6042; 1.187)) | inside | inside | inside (0.7: 0.1329 (0.5602; 1.319)) |
| 0.9 | 0.4452 (0.9136; 0.2592) → **0.2592** up | 0.6367 → **0.2592** up | 1.733 → **0.2592** up | 3.395 (0.8942; 0.3175) → **0.3175** up |
| 0.955 | 0.5442 (0.9136; 0.2592) → **0.2592** up | 0.8505 → **0.2592** up | 3.009 → **0.2592** up | 6.994 → **0.2592** up |

The class toy `default_inclusive_kernel(deuteron())` (ToyF2/ToyG1) binds at
two of the 32 points: x = 0.955, Q² = 2.5 and 5 (0.532 / 0.5285 → 0.4949,
|g₁|/F₁ = 0.835); x = 0.9 is inside there (0.5018 against 0.6525). On all
three deuteron kernels: f₁ f₂ g₁ g₂ b₃ b₄ Δ bit-identical at all 32 points; b₁,
b₂ bit-identical wherever nothing binds (24 / 30 / 24 points); post b₂ == 2x·b₁
exactly. ⁶Li and ⁷Li: all nine columns bit-identical on all three kernels.

### N2 — the A/B run list (`ab_runs.sh`, `ab_compare.py`; 300 events, seed 11)

| runs | npz arrays | npz meta | HepMC3 |
|---|---|---|---|
| 21 of P4's 22 CLI runs (`inc6_hul` exits 1 pre and post, A17), the 6 examples, the 24 inclusive ⁶Li plan × config runs | identical | identical (`b1_clipped_rate_fraction` absent) | identical |
| `incd` (`--isotope d --channel inclusive`, tensor-thirds, config 1, x_max 1.0) | **pre exit 1 (refused, P2) → post exit 0**: σ 495516.5512 pb, `b1_clipped_rate_fraction` 3.871134487837413e−7; banner "b1 spin-1 domain (…): binds on 3.871e-05 % of the rate, 59 cell(s) (59 upper, 0 lower; max unconfined b1/F1 = 6.993 at x = 0.955, Q2 = 1897)" | | |
| 24 inclusive d runs at `--x-max 0.85`: helicity-flip / apar (6) | identical | `b1_clipped_rate_fraction` (absent → 0.0), `rank2_input` (+ the confinement clause), `pol_sf_reach` (+ the g₁-bounds-b₁ clause); knob rows `b1_model`, `pol_sf` (reasons) | identical |
| … the 18 tensor-plan runs (azz, cos2phi, flip, tensor-thirds, transverse-tensor, tensor-flip) | identical | the above plus `pol_sf` ("not read under this run's fill" → `nnpdfpol`), `pol_sf_extrapolation` (label → `nearest`), `pol_sf_grid_q2_min` (NaN → 1.0), `pol_sf_below_grid_frac` (NaN → 0.0), `g1_clipped_rate_fraction` (NaN → 0.0); knob row `pol_sf` not-read → read | the run attributes `lipolgen.pol_sf` and `lipolgen.pol_sf_extrapolation` (label → `nnpdfpol` / `nearest`); nothing else |

The arrays of the deuteron runs are byte-identical because at x_max 0.85 the
window ends at the x = 0.794 row, below every binding cell (the key reads 0.0
there).

### N3 — builds at the default x_max 1.0 (`n3.py`; smallest φ-averaged density 1 + w_avg over the populated states and cells)

The alias names give the same rows as their plans. Pre (P2's tree): every
deuteron row refused (P2's M3 table). Post:

| d plan (CLI fill P_z 0.7, P_zz 0.6, P_e 0.7) | cfg 0 | cfg 1 | cfg 2 |
|---|---|---|---|
| helicity-flip / apar | builds; 0.2586 (apar+, m −1, x 0.871, Q² 660.1) | builds; 0.4549 (apar+, m −1, 0.871, 1897) | builds; 0.6974 (apar+, m −1, 0.871, 1897) |
| tensor-thirds / azz | builds; 0.8432 (azz+, m +1, 0.871, 660.1) | builds; 0.8398 (azz+, m +1, 0.871, 1897) | builds; 0.8400 (0.871, 1897) |
| transverse-tensor / cos2phi, tensor-flip / flip | builds; 0.8432 (m 0, 0.871, 660.1) | builds; 0.8398 | builds; 0.8400 |

**24 of 24 build.** ⁶Li: 24 of 24 build pre and post; every row of `n3.json`
(the CLI fill, the extreme fill and the pure states) identical pre and post;
the smallest densities are helicity-flip's 0.2117 / 0.1210 / 0.1965 (apar+,
m −1, x 0.955).

**The reach of the leading-order statement** (no action): (i) the extreme
fill P_z = P_zz = P_e = 1 (helicity-flip at `pzz_mode typed`, i.e. pure m = +1):
all 24 deuteron runs build; smallest densities helicity-flip 0.00462 / 0.2900 /
0.6362 (apar−, m +1, x 0.955 / 0.871 / 0.871), tensor-thirds and tensor-flip
0.8432 / 0.8398 / 0.8400, transverse-tensor 0.9957 / 0.9950 / 0.9950 (pre: 6 of
24 built, the transverse ones). (ii) The twelve pure states m ∈ {+1, 0, −1} ×
λ_e = ±1 at P_e = 1 × θ_S ∈ {0, 90°}, through the sampler's own check: none is
refused on any configuration; the smallest density is **0.00462** (config 0,
λ_e = +1, m = −1, θ_S = 0, x = 0.955, Q² = 733.6), 0.2900 and 0.6362 on configs
1 and 2. So at the most extreme states the window keeps a margin of 0.005 at
config 0's top row: the finite-γ g₂ term of the shipped `target_mass` kernel
and D(y) < 1 eat into the LO bound but do not cross it on the shipped grid.
⁶Li pure states: smallest 0.1272 / 0.0703 / 0.1778.

### N4 — `b1_clipped_rate_fraction` (`n45.py`; all 8 plan names × configs, x_max 1.0)

The definition (pipeline.hpp `B1RankTwoConfinement`, written before the
measurement): over the run's accepted cells, with w_c =
`Pipeline::cell_rate_weights_pb()`, the share of Σ w_c in the cells where, at
the cell centre, with the kernel's own F₁ and g₁ there (the sampler's cell
table) and the slot's unconfined b₁ (`InclusiveKernel::unconfined_b1`), F₁ > 0,
the window is not empty, and b₁ lies outside it.

| d cfg | fraction (identical for all 8 plan names) | cells (upper / lower) | rows | supervisor (joint.json) |
|---|---|---|---|---|
| 0 | 4.791951543386546e−7 | 40 (40 / 0) | 0.871, 0.955 | 4.791951543386554e−7, 40 |
| 1 | 3.871134487837413e−7 | 59 (59 / 0) | 0.871, 0.955 | 3.871134487837394e−7, 59 |
| 2 | 3.5361301668316917e−7 | 59 (59 / 0) | 0.871, 0.955 | 3.536130166831673e−7, 59 |

Reproduced: the same cells and halves; the shares agree to 1.7e−15 / 4.9e−15 /
5.1e−15 relative — the summation order (the library sums in cell order in
C++, the supervisor's script with numpy's pairwise sum). The largest
unconfined b₁/F₁: 5.422 (x 0.955, Q² 733.6) / 6.993 / 6.993 (x 0.955, Q² 1897).

### N5 — the domain after (`n7.py` post; full x_max 1.0 windows)

| d cfg | max(\|g₁\|/F₁ + b₁/(3F₁)) | excess over 1 | min b₁/F₁ | max b₁/F₁ |
|---|---|---|---|---|
| 0 / 1 / 2 | **1.0** / 1.0 / 1.0 (pre 2.7210 / 3.2446 / 3.2446) | 0 ULP | −0.02609 / −0.02994 / −0.03022 (unchanged) | 0.4713 / 0.4822 / 0.4822 (pre 5.422 / 6.993 / 6.993) |

P2's M4 table, post (config-1 kernel; b₁/F₁ at x = 0.5 / 0.7 / 0.8 / 0.85 /
0.9 / 0.955):

| Q² | post | (P2, pre) |
|---|---|---|
| 2.5 | −0.01219 / 0.03756 / 0.1026 / 0.1856 / **0.2592** / **0.2592** | 0.4452 / 0.5442 at 0.9 / 0.955 |
| 5 | −0.01418 / 0.04714 / 0.1355 / 0.2534 / **0.2592** / **0.2592** | 0.6367 / 0.8505 |
| 10.75 | −0.01624 / 0.05787 / 0.1743 / 0.3363 / **0.2592** / **0.2592** | 0.8816 / 1.278 |
| 52.39 | −0.01976 / 0.07891 / 0.2565 / **0.4399** / **0.2592** / **0.2592** | 0.521 / 1.465 / 2.43 |
| 150.6 | −0.02182 / 0.09267 / 0.314 / **0.5074** / **0.2592** / **0.2592** | 0.6559 / 1.916 / 3.419 |
| 167.3 | −0.02204 / 0.09411 / 0.3202 / **0.5129** / **0.2592** / **0.2592** | 0.6704 / 1.966 / 3.531 |
| 534.5 | −0.02459 / 0.1113 / 0.394 / **0.5526** / **0.3011** / **0.2592** | 0.8474 / 2.577 / 4.963 |
| 1897 | −0.02776 / 0.1329 / 0.4886 / **0.5666** / **0.3175** / **0.2592** | 1.078 / 3.395 / 6.994 |

(bold = confined). The cap 3(1 − 0.9136) = 0.2592 is where both nucleons'
g₁ are clipped to F₁ (A3); the confined b₁/F₁ is therefore no longer monotone
in x above x ≈ 0.85 — the shape of the domain, not of a model.

### N6 — HERMES (`m6_hermes.py`)

Every number of P2's M6 "post" column unchanged to every digit (`pre/m6.json`
== `post/m6.json`): `miller_shipped` χ²(b₁)/6 5.261863562 (4/6), A_zz 7.307981761
(3/6) / 5.549930644 (4/6) / 5.410395595; the kernel's own b₁ the same. At the
six bins nothing binds (the confined b₁ equals the unconfined one bit for bit):

| x | Q² | b₁/F₁ | \|g₁\|/F₁ | upper bound 3(1 − \|g₁\|/F₁) |
|---|---|---|---|---|
| 0.012 | 0.51 | 0.0046124 | 0.014457 | 2.9566 |
| 0.032 | 1.06 | 0.0059675 | 0.0080713 | 2.9758 |
| 0.063 | 1.65 | 0.0081705 | 0.022633 | 2.9321 |
| 0.128 | 2.33 | 0.0070219 | 0.090849 | 2.7275 |
| 0.248 | 3.11 | 0.0051999 | 0.20133 | 2.396 |
| 0.452 | 4.69 | −0.016422 | 0.35696 | 1.9291 |

(lower bound −1.5 everywhere).

### N7 — ⁶Li, record only (`n7.py`; nothing changed)

The default ⁶Li window's cells outside the JOINT half (the single-function
half holds everywhere, max b₁/F₁ 2.2904); config 0 has none:

| cfg | x | Q² | \|g₁\|/F₁ | b₁/F₁ | sum | rate share |
|---|---|---|---|---|---|---|
| 1 | 0.955 | 1536.1 | 0.28267 | 2.16735 | 1.00512 | 9.312e−12 |
| 1 | 0.955 | 1707.1 | 0.28267 | 2.22836 | 1.02545 | 7.834e−12 |
| 1 | 0.955 | 1897.2 | 0.28267 | 2.29043 | 1.04614 | 6.577e−12 |
| 2 | 0.955 | 1536.1 / 1707.1 / 1897.2 | 0.28267 | same | same | 1.044e−11 / 8.973e−12 / 7.699e−12 |

Totals 2.372e−11 (config 1) and 2.712e−11 (config 2), as the supervisor
found; no ⁶Li run is refused (N3: the smallest density 0.0703, a pure state
at P_e = 1).

### Labels, documents and pins (each with `J1 note (2026-11-12)`)

- New meta key `b1_clipped_rate_fraction` (inclusive deuteron runs only), the
  banner line, the `rank2_input` clause; `pol_sf` read on the inclusive
  deuteron (§1 N2).
- Tests re-taken: `tests/test_pipeline.cpp` "J1: …" (the slot equals
  `MillerB1` BEFORE the confinement, `unconfined_b1`; the confined value is the
  bound where it binds — at x 0.955, Q² 2.5 on the toy — and the factor
  3.0532 is checked on the unconfined slot); `python/tests/test_j1_deuteron_b1.py`
  (the same; the default run at x_max 1.0 again with its key pinned at
  3.871134487837413e−7, rel 1e−9 — a reproducibility tolerance chosen after
  the value was seen, said here; the pinned P2 refusal is now taken on the
  UNCONFINED deuteron kernel built through the public API — the shipped
  CT18NLO / NNPDFpol1.1 objects, `MillerB1`, the 1e−2 Δ, the option off —
  which reproduces P2's message: x 0.955, Q² 167.3, rank-2 −1.13, 24 cells);
  `python/tests/test_sf_backend.py`: the survivor of
  `test_the_neutron_tagged_pair_builds_since_a3_and_a_survivor_names_the_cell_and_the_cure`
  is that same unconfined deuteron kernel at helicity-flip config 1 (P2's
  pins −0.1021, −1.13, Q² 167.3 kept — chosen because it is the same physics
  the test pinned, now reached explicitly), and
  `test_the_cli_reports_it_as_a_configuration_error_not_a_traceback` takes
  its CLI refusal from the ⁶Li b₁ band ×2 row (`--b1-model cdks
  --b1-band-scale 2 --x-max 0.95`, refused at x = 0.871, the ⁶Li rank-2
  refusal the note says stands), cured at `--x-max 0.85` — the deuteron CLI
  run builds now; `python/tests/test_knob_provenance.py`: the
  `inclusive-d/tensor-thirds` spec back at x_max 1.0, `_spec_x_max`,
  `J1_XMAX_CELL`, the route-test hook and the `Spec.x_max` field removed (the
  file is the pre-P2 one bar two dated comments).
- New tests: C++ "J1 note: the spin-1 positivity confinement of the rank-2
  slot" (on/off, both halves, an unset and a set b₂, the empty window on a
  kernel whose eff_pol is 4.0, the non-spin-1 refusal, the factory on d and
  off on ⁶Li, the report on a d run and its absence on ⁶Li); Python
  `test_the_confinement_option_from_python`.
- Documents: `docs/USAGE.md` (the key defined next to g₁'s, "its absence
  means the confinement is not applied"; the refusal paragraph's J1 clause
  made true), `docs/PHYSICS_CHANNELS.md` (the rank-2 wiring row),
  `docs/OPEN_ITEMS_SOLUTIONS.md` §14, `docs/CONVENTIONS.md`, registry row 6's
  dated note (`run_2026-09-03/STATUS.md`), the CLI's `--x-max` help and refusal
  comment, comments in pipeline.hpp / sf.hpp; every edit line-neutral except
  the two appended blocks (`B1RankTwoConfinement` at the end of pipeline.hpp,
  `b1_rank2_confinement` at the end of pipeline.cpp) and the appended tests.

## 2. What was verified and how

- N1 and N2 are the A/B evidence (pre = `d86fbf2`'s build): ⁶Li, ⁷Li, tagged,
  coherent, the examples byte-identical in arrays, meta and HepMC3, the new
  key absent there; the deuteron's arrays at x_max 0.85 byte-identical; the
  deuteron's kernel moves only b₁, b₂ and only where the bound binds.
- The sampler's refusal (src/core/sampler.cpp `build_state`) is unchanged.
  Every sentence of its comment and message stays true of a deuteron run and
  of a kernel with the option on: the weight W is still never clamped ("NO
  CLAMP OF W" — what is confined is a structure function, to its own domain,
  inside the kernel, before any density is formed); "what this still catches
  is the RANK-2 term (the b1 band's x2 row, the toy b1 at the top cell …)"
  remains true for ⁶Li's band row, `--inclusive-b1` on the α tag and an
  unconfined caller kernel; the cure sentence names `--x-max` and the b₁
  model / band, both still the levers. No dated clause was needed.
- The two refusals still reachable on purpose are tested: the unconfined
  deuteron kernel (P2's numbers) and the ⁶Li band ×2 row.
- Line neutrality: `git diff --numstat` equal counts on every edited source
  and document but the appended blocks named above; the citation gates (§4).

## 3. Files

include/lipolgen/xsec.hpp, pipeline.hpp, sf.hpp; src/core/xsec.cpp,
pipeline.cpp; python/bindings.cpp; python/lipolgen/cli.py;
tests/test_pipeline.cpp; python/tests/test_j1_deuteron_b1.py,
test_sf_backend.py, test_knob_provenance.py; docs/USAGE.md,
PHYSICS_CHANNELS.md, OPEN_ITEMS_SOLUTIONS.md, CONVENTIONS.md,
open_items/run_2026-09-03/STATUS.md; this record and the P2b row of
open_items/run_2026-11-10/STATUS.md.

## 4. Gates

Foreground, one long job at a time, on the shared machine (load 1–6):

| gate | result | wall |
|---|---|---|
| build (`cmake --build build -j2`) | 0 warnings on the final tree (two intermediate `-Wmisleading-indentation` warnings on one-line `if`s fixed by braces; the full rebuild after the header change 135 s) | 25 s last incremental |
| `lipolgen_tests` (full) | **437 / 437 cases, 17 244 808 assertions**, 0 failed (P2: 436 / 17 244 738; +1 case, the option's) | 564 s |
| pytest, six parts (`pytest_part.sh`) | A 358 / 3 s; B 136 / 1 s; C 186; D 249 / 90 s; E 277 / 61 s; F 324 — **sum 1530 passed, 155 skipped, 0 failed** (P2: 1528 / 156: +1 the restored `x_max:1.0 -> 0.95` deuteron cell, +1 the new option test) | 203 / 116 / 81 / 411 / 292 / 189 s |
| bench subset (`LIPOLGEN_BENCH=1 … -k bench`) | **368 passed, 1 skipped** | 419 s |
| report `--all --plots --compare` (pre copy of REPORT.json) | **15 pass / 8 fail / 0 blocked; compare: 0 changes** (exit 0); the 19 regenerated figures differ from the committed ones only in the footer rows (pixel diff: the bottom 8 rows) and were restored with `git checkout -- docs/benchmarking/figures/*.png`; REPORT.json / REPORT.md differ only in `date`, `git_rev` and runtimes (left for the supervisor) | 109 s |
| report `--check` | OK | — |
| citations strict | 1243 / 95 / 7 / **0 broken** / 24 | 2.5 s |
| citations `--records-only` | 1432 citations in 95 records, **0 broken** | 0.3 s |
| citations `--audit-ranges` | 209, **0 REFUSED** (no range re-recorded) | 2.7 s |
| SPDX | 136 files, all carry the line | — |
| `test_doc_link_gate.py` | **112 passed** | 26 s |

## 5. Open points

1. **HEPMC3_CONVENTION clause for P4 (J6).** The new key is NOT a HepMC3 run
   attribute (g₁'s share is not one either), so no clause is required; but
   the convention's sentence on `pol_sf` should know that a deuteron file now
   carries `lipolgen.pol_sf` = the backend name on the tensor plans too. If P4
   wants it, the sentence to add after the list of keys: "On an inclusive
   deuteron run `pol_sf` (and `pol_sf_extrapolation`) names the backend under
   every fill, since g₁ bounds the deuteron's rank-2 slot (J1 note,
   2026-11-12); `meta["b1_clipped_rate_fraction"]`, the share of the rate that
   bound binds on, is an npz key only, like `g1_clipped_rate_fraction`."
2. **`--pol-sf` is read on every inclusive deuteron run** — not foreseen by
   the brief ("the meta differ by exactly the new key and the rank2_input
   clause"): the tree's own knob gate (`test_knob_provenance.py`, "moved →
   read") failed on `pol_sf:nnpdfpol -> toy` on the deuteron's tensor-thirds
   spec, because g₁ now decides b₁ on the binding cells. `pol_sf_is_read`
   says so for the inclusive deuteron (config-level: also at an x_max whose
   window has no binding cell, where g₁ is still evaluated against b₁ on every
   cell but moves nothing — the 0.85 runs of N2); its meta and HepMC3
   consequences are listed in N2. The supervisor may prefer a run-level
   reading ("read where it binds"); it would cost a Pipeline argument to a
   config-level function.
3. ⁶Li's three joint-half cells (N7) stay as recorded — the modelling
   question the note files for the author (reading the toy as a
   scale-independent b₁/F₁) would remove them on both ions.
4. The class toy deuteron kernel (`default_inclusive_kernel(deuteron())`) is
   confined too (two of P2's 32 points move); no reference JSON pins it (the
   C++ suite and the port gates pass).
5. The confined deuteron b₁/F₁ is not monotone above x ≈ 0.85 (N5) — the cap
   follows A3's clipped g₁. A reader of a deuteron A_zz at x > 0.85 should
   read it as "the domain's edge", not as a model prediction.
