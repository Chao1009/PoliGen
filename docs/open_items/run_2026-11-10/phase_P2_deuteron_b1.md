<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-11-10 — phase P2: the deuteron's own b₁ slot (decision J1)

Worker P2, 2026-11-10, on `claude/run-2026-10-26` at `1c98f62` (phase P1 of
this run committed; this phase uncommitted, the supervising session commits).
Decision: `DECISIONS.md` J1 — **the AUTHOR's decision** ("1. go with your
recommendation", on the recommendation of `../run_2026-11-03/SUMMARY.md` item
1). J1 has two halves: the deuteron's slot (implemented here) and the
re-opening of the toy b₁'s high-x behaviour (MEASURED here, §1 M3–M5; no
treatment is chosen or implemented — no clip, no taper change, no x_max
change; the supervisor decides by a dated note under J1). Every number below
was measured in this phase with `source env.sh` on this machine; scripts,
logs and the pre/post outputs are in the session scratch `run_1110/p2/`
(`scripts/m1_dump.py`, `m3_refusals.py`, `m45.py`, `m6_hermes.py`,
`ab_runs.sh`, `ab_compare.py`; `pre/`, `post/`, `ab/pre`, `ab/post`,
`ab_compare.txt`, `ab_summary.txt`, `ab_sigma.txt`, `ab_events.txt`, `py/`,
the gate logs). **pre** = `1c98f62`'s build, every pre measurement taken
before the first edit; **post** = this tree. *Supervisor's addition,
2026-11-12:* the worker was stopped by the usage limit during its last gate
(2026-11-10), and a reboot of the machine on 2026-11-12 wiped the session
scratch — the pre/post outputs and logs named above no longer exist; the
scripts were restored from the session transcript. Every number below is
the worker's, measured before that; §4 says which gate ran when.

## 0. What changed, in one paragraph

`default_inclusive_kernel` (src/core/pipeline.cpp) keeps its one spin test
`kernel_fills_rank2` and now keys WHICH camp fills the spin-1 slot on what the
ion IS (`spin1_rank2_camp`: A = 2, Z = 1 → the deuteron; A = 6, Z = 3 → ⁶Li;
any other spin-1 ion THROWS, naming the missing mapping). The deuteron's b₁
slot is `MillerB1` itself, per nucleon, unmapped — a `shared_ptr` captured by
value in the closure, as the Cdks branch does (no raw-`this` closure, nothing
dangles); ⁶Li keeps the function-local `static Li6B1(MillerB1)` and its
explicit-transfer branch line for line; the Δ slot is the same for both. On
the deuteron a non-Miller `model` and a moved `b1_transfer` are REFUSED with a
sentence that says why (the refused value printed `%.17g`); the default is
accepted and the meta key `li6_b1_rank2_transfer` carries the label `not read
on isotope d -- the deuteron's rank-2 slot is MillerB1 itself, per nucleon,
unmapped (J1, 2026-11-10)` (A17's rule; `InclusivePolarizationFactors::
li6_b1_rank2_transfer_label` / `meta_labels()`, written by `run_provenance`
and by the npz `meta`); the banner's "inclusive polarization factors" line and
`meta["rank2_input"]` say the same on a deuteron run. `validate()` is
unchanged (`--b1-model cdks` / `li6-convolution` stay refused off ⁶Li). The
measured consequence: every inclusive deuteron rank-2 number ×3.0531989383
(= 1/0.32753), ⁶Li / ⁷Li / tagged / coherent byte-identical, and **every
default inclusive deuteron run — every CLI plan, configs 0, 1, 2 — is now
refused at x_max 1.0** by the rank-2 term on the whole x = 0.955 row (§1 M3),
a share 1.6 × 10⁻⁹ of the run's rate (M5).

## 1. Every number or label that moved, before → after

### M1 — kernel tables (`scripts/m1_dump.py`; `pre/m1.json`, `post/m1.json`, `post/m1_check.txt`)

x ∈ {0.012, 0.05, 0.1, 0.3, 0.5, 0.7, 0.9, 0.955} × Q² ∈ {2.5, 5, 100, 1897}
(all 32 points tabulated: the kernel evaluates anywhere; the in-window subset
is a subset), nine columns f₁ f₂ g₁ g₂ b₁ b₂ b₃ b₄ Δ, on three kernels per
ion: `shipped_inclusive_kernel(ion)` (D7 grids), the class toy
`default_inclusive_kernel(ion)`, and the `Pipeline`'s own kernel.

| ion | kernel | b₁, b₂ post/pre | f₁ f₂ g₁ g₂ b₃ b₄ Δ |
|---|---|---|---|
| d | shipped / toy / pipeline | **3.053198938300956** at all 32 points (min = max), = 1/((2/6) × 0.982576); post b₁ == `MillerB1().b1` bit for bit | bit-identical (ratio exactly 1) |
| ⁶Li | shipped / toy / pipeline | bit-identical | bit-identical |
| ⁷Li | shipped / toy / pipeline | bit-identical | bit-identical |

### M2 — the A/B run list (`scripts/ab_runs.sh`, `ab_compare.py`; `ab_summary.txt`, `ab_sigma.txt`, `ab_events.txt`)

Run 2026-11-03 P4's list (22 CLI runs + 6 examples; 300 events, seed 11)
plus every CLI plan name (`apar azz cos2phi flip helicity-flip tensor-flip
tensor-thirds transverse-tensor`; the first four are aliases of the last
four, measured identical) × configs 0, 1, 2 on the inclusive deuteron at
`--x-max 0.85` and on inclusive ⁶Li at the default 1.0 (24 + 24 runs).

| runs | npz arrays + meta | HepMC3 |
|---|---|---|
| 21 of P4's 22 CLI runs (`inc6_hul` exits 1 pre and post, the A17 refusal) and the 6 examples | byte-identical | byte-identical |
| `incd` (`--isotope d --channel inclusive`, default plan tensor-thirds, config 1, x_max 1.0) | **pre exit 0 → post exit 1**: refused, m = 1 at x = 0.955, Q² = 167.3, 1 + w_avg = −0.1302, rank-2 −1.13, 24 cells | — |
| 24 inclusive ⁶Li plan × config runs | byte-identical | byte-identical |
| 24 inclusive d plan × config runs at x_max 0.85 | arrays that move: the event kinematics of the events whose cell draw moved (`x q2 y w2 nu e_prime theta_e eta_e kp cell`, plus `phi`, `struck_pdg` on helicity-flip c0) — 17 to 80 of 300 events; `category` and `m` identical; meta keys that move: `sigma_pb`, `sigma_gen_mb`, `sigma_per_category_pb`, `rank2_input` (the J1 clause appended), `li6_b1_rank2_transfer` (0.982576 → the `not read …` label), and the `b1_model` knob row's reason (it carries `rank2_input`); nothing else | the run attribute `lipolgen.li6_b1_rank2_transfer` (same label), `GenCrossSection`, the moved events |

The cross sections (`ab_sigma.txt`; the spin-blind total cannot move on a
fill with zero net alignment, the per-category ones move by the tensor term):
tensor-thirds σ unchanged to every digit, per category ×0.999420 / 0.999420 /
1.001160 (c0), ×0.999595 / 0.999595 / 1.000809 (c1), ×0.999731 / 0.999731 /
1.000538 (c2); helicity-flip σ ×0.99960394 / 0.99972370 / 0.99981626 (c0 /
c1 / c2); transverse-tensor ×1.00029012 / 1.00020242 / 1.00013462;
tensor-flip ×0.99985491 / 0.99989878 / 0.99993268.

### M3 — the exposure at the default x_max 1.0 (`scripts/m3_refusals.py`; `pre/m3.json`, `post/m3.json`, `post/m3_hf_c0.txt`)

The numbers are the sampler's refusal message (the FIRST refused pure spin
state m and cell in build order; "cells" = refused cells of that state); the
last column is the largest `--x-max` on a 0.01 step that builds. P_z 0.7,
P_zz 0.6, P_e 0.7 (the CLI defaults). The alias plans give the same rows as
their names (`apar` = helicity-flip, `azz` = tensor-thirds, `cos2phi` =
transverse-tensor, `flip` = tensor-flip).

| d plan | cfg | pre | post: m, Q², 1 + w_avg (vector, rank-2), cells | post x_max |
|---|---|---|---|---|
| helicity-flip | 0 | refused: m −1, Q² 534.5, −0.05837 (−0.5213, −0.537), 4; x_max 0.95 | m 1, 186, −0.005743 (+0.1664, −1.172), 14 | **0.87** |
| helicity-flip | 1 | refused: m −1, 1536, −0.02078 (−0.3019, −0.7189), 3; 0.95 | m 1, 167.3, −0.1021 (+0.02803, −1.13), 24 | 0.95 |
| helicity-flip | 2 | builds | m 1, 167.3, −0.1189 (+0.01118, −1.13), 24 | 0.95 |
| tensor-thirds | 0 / 1 / 2 | builds | m 1, 167.3, −0.1314 / −0.1302 / −0.1301 (0, −1.131 / −1.13 / −1.13), 15 / 24 / 24 | 0.95 |
| transverse-tensor | 0 / 1 / 2 | builds | m 0 (θ_S = 90°), same numbers as tensor-thirds | 0.95 |
| tensor-flip | 0 / 1 / 2 | builds | m 0 (θ_S = 90°), same numbers as tensor-thirds | 0.95 |

So 4 of 24 deuteron runs were refused pre, **24 of 24 post**; the lowest
refused x is 0.955 on every message. The cell grid is FIXED (rows
x = 0.794, 0.871, 0.955; `--x-max` drops rows, it does not move them), so
x_max 0.95 = "without the 0.955 row" and 0.87 = "without the 0.871 row" too.
helicity-flip config 0 needs 0.87 because ANOTHER state is refused one row
lower: at x_max 0.95 the message is m = −1 at x = 0.871, Q² = 660.1,
1 + w_avg = −0.02616 (vector −0.5846 at A₁ 0.8429, rank-2 −0.4416) — the
vector and rank-2 terms together, each alone above −1.
**⁶Li: all 24 plan × config runs build at 1.0, pre and post; the 24 rows of
`pre/m3.json` and `post/m3.json` are identical.**

### M4 — the positivity map of the toy (post; `scripts/m45.py`, `post/m45.json`)

**The condition, from `InclusiveKernel::tensor_amplitudes`** (the default
`tensor_gamma` = false path) and `phi_averaged_density`: for a pure state m of
the J = 1 ion at spin angle θ_S, with an unpolarised lepton (λ_e P_e = 0, so
the vector term is zero), the φ-averaged density is 1 + w with

w(m) = TENSOR_LL_SIGN · Q_NN(m) · P₂(cos θ_S) · [b₁ + (1−y)/(x y²) b₂] / [F₁ + (1−y)/(x y²) F₂],

TENSOR_LL_SIGN = −1, Q_NN(m) = (3m² − 2)/3 (+1/3 at m = ±1, −2/3 at m = 0).
With the kernel's b₂ = 2x b₁ and its own F₁ = F₂/(2x(1+R)) this is
**w(m) = −Q_NN(m) P₂ (b₁/F₁) ρ**, ρ(y, R) = [1 + 2(1−y)/y²] / [1 + 2(1+R)(1−y)/y²]
∈ (1/(1+R), 1], ρ → 1 as y → 1. Fully aligned (θ_S = 0): m = ±1 needs
**b₁/F₁ ≤ 3/ρ**, m = 0 needs **b₁/F₁ ≥ −3/(2ρ)**; for every y (sup ρ = 1):
**−3/2 ≤ b₁/F₁ ≤ 3** — the structure-function-level bound, the rank-2
analogue of |g₁| ≤ F₁. (At θ_S = 90°, the transverse plans: m = 0 carries
−(1/3)ρ b₁/F₁, the same upper bound.) The Δ term enters a₂ (cos 2φ) only,
not the φ-averaged density. Check at the first refused cell (x 0.955, Q² 167.3):
b₁/F₁ = 3.531, w(+1) = −1.1302 → ρ = 0.960, the bound 3/ρ = 3.125.

**On the shipped kernels** (config-1 grid; "first negative" = the lowest-x
cell of that Q² row where a pure aligned state's density is ≤ 0, its m and
rank-2 weight; Q² = 2.5 / 5 are kernel scans of the y-free bound on
x = 0.01 … 0.99, they have no cells at high x):

| ion | Q² row | top cell x | first negative | b₁/F₁ at x = 0.5 / 0.7 / 0.8 / 0.85 / 0.9 / 0.955 |
|---|---|---|---|---|
| d | 2.5 | — | none (bound holds to 0.99) | −0.01219 / 0.03756 / 0.1026 / 0.1856 / 0.4452 / 0.5442 |
| d | 5 | — | none | −0.01418 / 0.04714 / 0.1355 / 0.2534 / 0.6367 / 0.8505 |
| d | 10.75 | 0.550 | none | −0.01624 / 0.05787 / 0.1743 / 0.3363 / 0.8816 / 1.278 |
| d | 52.39 | 0.871 | none | −0.01976 / 0.07891 / 0.2565 / 0.521 / 1.465 / 2.43 |
| d | 150.6 | 0.871 | none | −0.02182 / 0.09267 / 0.314 / 0.6559 / 1.916 / 3.419 |
| d | 167.3 | 0.955 | **x 0.955, m = ±1, w −1.130** | −0.02204 / 0.09411 / 0.3202 / 0.6704 / 1.966 / 3.531 |
| d | 534.5 | 0.955 | **x 0.955, w −1.630** | −0.02459 / 0.1113 / 0.394 / 0.8474 / 2.577 / 4.963 |
| d | 1897 | 0.955 | **x 0.955, w −2.322** | −0.02776 / 0.1329 / 0.4886 / 1.078 / 3.395 / 6.994 |
| ⁶Li | 2.5 | — | none | −0.003994 / 0.0123 / 0.03362 / 0.06078 / 0.1458 / 0.1782 |
| ⁶Li | 167.3 | 0.955 | none | −0.007218 / 0.03082 / 0.1049 / 0.2196 / 0.6438 / 1.157 |
| ⁶Li | 1897 | 0.955 | none (the window's maximum, 2.291 < 3) | −0.009091 / 0.04352 / 0.16 / 0.3531 / 1.112 / 2.291 |

Every one of the 68 Q² rows of the config-1 grid is in `post/m45.json`: on
the deuteron the pure-state density is negative on EVERY cell of the x = 0.955
row (every row from Q² = 167.3 up; that row has no cell below) and on no cell
of any other row; on ⁶Li on none. b₁/F₁ rises with Q² at fixed x because
`MillerB1` has no Q² dependence while CT18NLO's F₁ falls at high x as Q² grows
(x = 0.955: 0.544 at Q² 2.5 → 6.99 at 1897 on the deuteron).

**What the toy is up there.** `MillerB1` = 0.5 × Miller's per-deuteron table
(`B1_MILLER_TABLE_TO_PER_NUCLEON`), digitized on x = 0.010 … 0.900;
`toy_b1` reads it with `DigitizedTable::interp_tapered` (x_end 1, power 3):
constant below the table (0.0646795 per nucleon at x = 0.001 and 0.005), and
above x = 0.900 the table's last value times ((1 − x)/0.1)³ — 1.5622e−4 at
0.900, 8.0004e−5 at 0.92, 1.4239e−5 at 0.955 (= 1.5622e−4 × 0.0911). So the
whole x = 0.955 row reads this code's own taper, not Miller's curve, and the
taper (∝ (1−x)³) falls more slowly than CT18NLO's F₁ there at large Q².

### M5 — the rate at stake (post; `post/m45.json` "M5")

Share of the spin-blind cell rate (`InclusiveSampler::cell_xsec_pb` of the
full x_max 1.0 window, which equals `Pipeline.cell_rate_weights_pb` on the
inclusive channel — checked equal on a built run) in (i) the cells some plan
state refuses at x_max 1.0 (my own per-cell evaluation over every populated
(category, m) state of the plan with `amplitudes_at`; it reproduces the
messages' per-state counts, e.g. 14 for helicity-flip c0 m = +1, 24 for
tensor-thirds c1), and (ii) the cells where a pure aligned state's density is
negative (θ_S = 0 and θ_S = 90° give the same cells):

| ion, cfg | window cells | (i) refused cells, share — tensor-thirds / transverse / tensor-flip | (i) helicity-flip | removed by the largest building x_max | (ii) pure-state-negative cells, share | g₁ clip share (A3, for scale) |
|---|---|---|---|---|---|---|
| d, 0 | 1950 | 15 (row 0.955), 1.677e−9 | 16 (+1 at x 0.871), 3.208e−9 | 1.677e−9 (0.95); helicity-flip 4.792e−7 (0.87) | 15, 1.677e−9 | 3.5–3.9e−7 (P2 of run 11-03) |
| d, 1 | 3053 | 24, 1.658e−9 | 24, 1.658e−9 | 1.658e−9 | 24, 1.658e−9 | |
| d, 2 | 3498 | 24, 1.563e−9 | 24, 1.563e−9 | 1.563e−9 | 24, 1.563e−9 | |
| ⁶Li, 0 / 1 / 2 | 1947 / 3051 / 3498 | 0 | 0 | 0 (builds at 1.0) | 0 | |

The whole x = 0.955 row carries 1.6 × 10⁻⁹ of the rate (on ⁶Li too: 1.676e−9 /
1.655e−9 / 1.563e−9); x > 0.83 (the 0.871 and 0.955 rows) carries 4.79e−7 /
3.87e−7 / 3.54e−7. Any treatment of the high-x toy touches ~10⁻⁹ of the
deuteron rate (and the 0.871 row's 4.8e−7 only on helicity-flip config 0).

### M6 — HERMES (`scripts/m6_hermes.py`; `pre/m6.json`, `post/m6.json`)

| quantity | pre | post |
|---|---|---|
| `miller_shipped` (the headline, `MillerB1()`): χ²(b₁)/6, bins | 5.261863562, 4/6 | 5.261863562, 4/6 (unchanged) |
| … A_zz on the shipped F₁ᵈ / on MSTW-R1990 / with R1998 | 7.307981761 (3/6) / 5.549930644 (4/6) / 5.410395595 | unchanged |
| the shipped deuteron kernel's OWN b₁ (`shipped_inclusive_kernel(deuteron()).tables(x, Q²).b1`) / `MillerB1` per bin | 0.327525333 at all six bins | **1** (equal bin for bin) |
| … χ²(b₁)/6, bins | 14.557823106, 2/6 | **5.261863562, 4/6** |
| … A_zz on the shipped F₁ᵈ / on MSTW-R1990 | 15.678316250 (2/6) / 14.681725308 (2/6) | **7.307981761 (3/6) / 5.549930644 (4/6)** |
| row status | FAIL | FAIL |

After J1 the kernel's own b₁ IS `MillerB1`, so the headline and what a
default deuteron run carries coincide, to every digit. The harness never read
the kernel's b₁, so the report row does not move (§2, the `--compare`).

### Labels, documents and pins re-taken (each with `J1 (2026-11-10)`)

- `meta["li6_b1_rank2_transfer"]` / HepMC3 `lipolgen.li6_b1_rank2_transfer`
  on an inclusive deuteron run: `0.982576` → the `not read on isotope d …`
  label; banner: "; 6Li rank-2 transfer 0.982576 (…), which the deuteron's
  class-level b1 reads too" → "; rank-2 slot MillerB1 per nucleon, unmapped:
  the 6Li rank-2 transfer is not read (J1, 2026-11-10; 0.32753 x MillerB1
  before)"; `meta["rank2_input"]` (and the `b1_model` knob row's reason) on a
  deuteron run gains "; on the deuteron b1_func is MillerB1 itself, per
  nucleon, unmapped, and the 6Li rank-2 transfer is not read (J1, 2026-11-10)".
- `python/tests/test_sf_backend.py::test_the_neutron_tagged_pair_builds_since_a3_and_a_survivor_names_the_cell_and_the_cure`:
  the inclusive-deuteron refusal message re-taken, `1 + w_avg = -0.02078` →
  `-0.1021`, `-0.7189 from the RANK-2 term` → `-1.13 …` (+ `Q2 = 167.3`).
- `python/tests/test_knob_provenance.py`: the `inclusive-d/tensor-thirds`
  spec's BASE is refused at 1.0 since J1, so 43 of its cells skipped and
  three tests failed; it runs at x_max 0.95 again through `Spec(x_max=…)` and
  `_spec_x_max` (the A1-era helper P1 of run 11-03 had and P2 removed, with a
  dated docstring), the route test likewise for that one case, and the one
  cell that cannot exist then, `x_max:1.0 -> 0.95` on that spec, is skipped
  by name with its reason (`J1_XMAX_CELL`). No check loosened: every other
  deuteron cell runs as before J1.
- New: `tests/test_pipeline.cpp` "J1: the inclusive deuteron's rank-2 slot is
  MillerB1 itself; 6Li keeps Li6B1(MillerB1)" (bitwise b₁/b₂/Δ, the factor,
  the refusals, the old reading built explicitly, the throw on another spin-1
  ion, meta / banner / `run_provenance` on d and on ⁶Li);
  `python/tests/test_j1_deuteron_b1.py` (6 tests: the same from Python, the
  shipped kernel = the pipeline's, the meta label, and M3's tensor-thirds
  refusal pinned).
- No test whose PURPOSE is the old mapping exists (searched: every
  `default_inclusive_kernel(deuteron())` caller reads f₂ / g₁ / validate
  paths); none needed `Li6B1` passed explicitly. No C++ pin moved.
- Documents (dated clauses, line-neutral): `docs/PHYSICS_CHANNELS.md` (the
  "for ANY spin-1 ion" row), `docs/OPEN_ITEMS_SOLUTIONS.md` §14 "Adjacent
  defects" (the bullet CLOSED by J1), `docs/USAGE.md` (the refusal paragraph
  of the top x cell), `docs/CONVENTIONS.md` (the transfer is read by ⁶Li
  alone), `validation/benchmarks/README.md` row 3 and convention 8 (the
  0.32753 / 14.56 clauses become history), the harness's J1 NOTE (after the
  A5 NOTE; its first 160 lines byte-identical), registry row 6 of
  `run_2026-09-03/STATUS.md` (a dated note; no label moves), comments in
  `pipeline.cpp` / `pipeline.hpp` / `constants.hpp` / `sf.hpp`, the CLI's
  `--x-max` help and its refusal comment, `shipped_inclusive_kernel`'s
  docstring.

## 2. What was verified and how

- M1 and M2 above are the A/B evidence: ⁶Li, ⁷Li, tagged, coherent and the
  six examples byte-identical in arrays, meta and HepMC3; the deuteron's
  b₁/b₂ moved by exactly one factor and nothing else in its tables.
- The old deuteron reading stays constructible: from C++ as
  `InclusiveKernel(DEUTERON(), Options{b1_func = Li6B1(MillerB1)})` (the new
  C++ test), and from Python as `InclusiveKernel(deuteron(), o)` with
  `o.b1_func = lambda x, q2, f1: Li6B1(MillerB1()).b1(x, q2, f1)` (the new
  pytest). No flag was added.
- Lifetime: the deuteron closure owns its `MillerB1` by value; the ⁶Li
  default still uses the function-local static the comment there describes.
- Line neutrality: every edited source and document has equal added and
  removed line counts (`git diff --numstat`), bar the appended tests, the new
  files and the harness's NOTE (no citation reaches past its line 17); no
  range fingerprint needed re-recording (the citation gates, §4).
- Report: `make_report.py --all --plots --compare REPORT_pre.json` — see §4.

## 3. Files

src/core/pipeline.cpp; include/lipolgen/pipeline.hpp, constants.hpp, sf.hpp;
python/bindings.cpp; python/lipolgen/__init__.py, cli.py;
tests/test_pipeline.cpp; python/tests/test_j1_deuteron_b1.py (new),
test_knob_provenance.py, test_sf_backend.py;
validation/benchmarks/t2_hermes_b1_table2.py, README.md;
docs/PHYSICS_CHANNELS.md, OPEN_ITEMS_SOLUTIONS.md, USAGE.md, CONVENTIONS.md,
open_items/run_2026-09-03/STATUS.md; this record and the P2 row of
open_items/run_2026-11-10/STATUS.md.

## 4. Gates

| gate | result | wall | run by |
|---|---|---|---|
| build | up to date with the tree (the worker's last `cmake --build build -j2`, 2026-11-10); its warning count is not recoverable (log lost) — phase P2b rebuilds the same files and reports it | — | worker |
| `./build/lipolgen_tests` | 436 cases, 17 244 738 assertions, all passed (435 / 17 244 595 before: + the J1 case) | 419 s | worker, final tree |
| pytest part A / B | 358 passed, 3 skipped / 136 passed, 1 skipped | 205 s / 131 s | worker, before the re-pins of `test_knob_provenance.py` and `test_sf_backend.py` (parts D–F); no file of A or B, and no source, changed after |
| pytest part C / D / F | 186 passed / 248 passed, 91 skipped / 323 passed | 71 s / 351 s / 162 s | worker, final tree |
| pytest part E | 277 passed, 61 skipped | 429 s | supervisor, 2026-11-12 (the worker's run finished unread) |
| pytest, sum | **1528 passed, 156 skipped, 0 failed** (1523 / 155 before: + 6 tests of `test_j1_deuteron_b1.py`, one knob cell passed → skipped by name) | | |
| bench subset | 368 passed, 1 skipped | 527 s | supervisor, 2026-11-12, on a loaded machine (a first attempt was killed by its 590 s limit at 366 of 369 and is not counted) |
| report `--all --plots --compare` (baseline: the committed `REPORT.json`) | **0 changes**; 15 pass / 8 fail / 0 blocked; `--check` OK; the regenerated files differed by stamps only and were restored | 145 s | supervisor, 2026-11-12 |
| citations | strict 1243 references / 95 ranges / 7 external / 0 broken / 24 allow-listed; `--records-only` 1432 citations in 94 records, 0 broken; `--audit-ranges` 209, 0 REFUSED | | supervisor, 2026-11-12 |
| SPDX; doc-link gate; `git diff --check` | 136 tracked files + the two new ones carry the line; 112 passed; clean | 39 s | supervisor, 2026-11-12 |

## 5. Open points

*Supervisor's addition, 2026-11-12:* item 1 is decided by the dated note
under J1 in `DECISIONS.md` — the FULL spin-1 domain (b₁ yields to g₁ too),
on the deuteron's slot, implemented in phase P2b; the worker's
recommendation below is kept as written.

1. **The exposure, for the supervisor's high-x decision under J1.** Facts:
   every default inclusive deuteron run (24 of 24 plan × config) is refused at
   x_max 1.0 (4 of 24 before J1); the refused cells are exactly the x = 0.955
   row, 15 / 24 / 24 cells, 1.6 × 10⁻⁹ of the rate; the pure aligned states
   go negative there and nowhere else (b₁/F₁ up to 6.99 against the bound 3);
   that row reads this code's (1−x)³ taper above Miller's last digitized
   point x = 0.900 on a Q²-independent table, against a Q²-evolving CT18NLO
   F₁; ⁶Li stays inside the bound everywhere (max 2.29) and builds; x_max
   0.95 builds every deuteron run but helicity-flip config 0 (0.87, a joint
   vector + rank-2 cell at x = 0.871, each term alone above −1).
   **Recommendation (not implemented): (b), a visible positivity bound on the
   rank-2 term, A3's treatment: −3/2 ≤ b₁/F₁ ≤ 3 per nucleon where the
   kernel's b₁ is evaluated, the clipped rate share in `meta` and on the
   banner as g₁'s is.** Reason: it is the exact, y- and convention-free
   bound of §1 M4 (the rank-2 counterpart of |g₁| ≤ F₁, which A3 already
   chose for the same kind of defect); it binds only on cells that are the
   toy's own taper (x = 0.955 on the deuteron, 1.6e−9 of the rate) and on none
   of ⁶Li's, so ⁶Li stays bit-identical and every deuteron number below the
   top row too; option (a) leaves `--isotope d` refused out of the box on
   every plan and makes a smaller x_max the de-facto default, which A4 said
   only hides a violation; a taper or x_max change would be a choice among
   many unmeasured shapes. What (b) does NOT lift: helicity-flip config 0's
   x = 0.871 cell (−0.5846 vector, −0.4416 rank-2: a JOINT violation at
   P_e P_z = 0.49, each term inside its own bound); that one would stay a
   refusal under A4 (cured by `--x-max 0.87`), or need the joint bound
   1 + w₂(m) ≥ |w₁(m)| at |λ_e| = 1 — a larger design, not recommended now.
2. **`--b1-model cdks` on the deuteron** would be well defined after J1
   (`CdksB1` unmapped in the deuteron's own slot); it stays refused, in
   `validate()` and now also in `default_inclusive_kernel` (an API call with
   a non-Miller model on the deuteron throws instead of silently building the
   ⁶Li mapping). Enabling it is a separate decision.
3. `docs/HEPMC3_CONVENTION.md` says the three A1 keys are written "where the
   inclusive run reads them"; on a deuteron run `li6_b1_rank2_transfer` is
   now written as a `not read …` label (the same file allows labels:
   "labels included"). The document is not edited (rule 3); J6 / phase P4
   may want one clause.
4. The USAGE refusal table of the top x cell is D7's (2026-10-29, marked so
   by V-CODE-5); the J1 clause beside it states today's deuteron rows; a full
   re-measured table is left to whoever acts on item 1.
