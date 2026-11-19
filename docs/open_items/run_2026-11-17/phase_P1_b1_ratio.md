<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-11-17 — phase P1: the digitized b₁ read as the ratio b₁/F₁ at a reference scale (K1)

Worker P1 (Opus, under the supervising session), on `claude/run-2026-10-26`
at base `a1f2b4d` (tree clean at start; this phase uncommitted, the
supervisor commits). Started 2026-11-17, stopped by the usage limit, resumed
and finished 2026-11-19 (the date of every measurement below unless said).
What it implements: **K1** of `DECISIONS.md` — **the author's decision**
("go with your recommendation"): the digitized deuteron b₁ tables are read
as b₁(x, Q²) = b₁ᵗᵃᵇˡᵉ(x) · F₁(x, Q²)/F₁(x, Q₀²); its details (a)–(g)
(Q₀² = 2.5 GeV², the kernel's own F₁, the scope, the run surface only, one
knob, b₂ in step, F₁(x, Q₀²) ≤ 0 → 0) are **the supervisor's choices** under
it, implemented as written. **pre** = the build of `a1f2b4d`, snapshotted
before the first edit (A/B set, `REPORT.json` copy, a copy of the Python
package and libraries); **post** = this tree. Scripts and raw logs are in the
session scratch `run_1117/p1/` (volatile; every number a reader needs is
here). Dated comments of this phase: `K1 (2026-11-17)`.

## 0. What changed, in one paragraph

One kernel option, `InclusiveKernel::Options::b1_reference_q2` (xsec.hpp;
**default 0 = the frozen table**, so a class-level kernel is bit for bit the
pre-K1 one), applied in `InclusiveKernel::tables` on the spin-1 slot only,
after the slot is filled and **before** the J1 note's confinement: b₁ ←
b₁ · F₁(x, Q²)/F₁(x, Q₀²) with the kernel's own per-nucleon F₁ at both scales
(`b1_reference_factor`; 0 where F₁(x, Q₀²) ≤ 0), a set b₂ by the same factor,
an unset b₂ recomputed as 2x·b₁; `unconfined_b1` returns the scaled slot, so
`b1_rank2_confinement` reports the share of what `tables()` actually clips.
The constructor refuses a negative / non-finite value and a positive one off
spin 1. One constant `B1_REFERENCE_Q2 = 2.5` (constants.hpp). The run
surface: `PipelineConfig::b1_reference_q2 = B1_REFERENCE_Q2`, CLI
`--b1-ref-q2` (Python `make_config(b1_reference_q2=…)`), read where a
digitized table is the rank-2 input (inclusive 6Li / d at `--b1-model miller`
or `cdks`, and tagged-6Li-alpha with `--inclusive-b1` through a new
`StruckClusterOptions::b1_reference_q2`), labelled not-read elsewhere and a
moved value refused there (A17); with an LHAPDF unpolarised grid a value
outside the grid's own Q² range is refused (new virtual
`UnpolSF::grid_q2_range`, answered by `LhapdfSF`). `default_inclusive_kernel`
gains a trailing `b1_reference_q2 = 0` (refused positive on
`li6-convolution`); `lipolgen.shipped_inclusive_kernel` gains the keyword,
default 0 (decision (d): factories keep the frozen table unless asked).
`meta["b1_reference_q2"]` (npz only), a banner line, a `rank2_input` clause,
a knob row. The HERMES row's headline is re-pointed to the run-surface
kernel's b₁ (convention 8); the frozen table is a recorded sub-row. Measured:
b₁/F₁ is Q²-independent at fixed x to 2.2e−16 under K1; the deuteron's
confinement share falls from 4.79e−7 / 3.87e−7 / 3.54e−7 (40 / 59 / 59 cells)
to 1.68e−9 / 1.66e−9 / 1.56e−9 (15 / 24 / 24 cells, x = 0.955 only); ⁶Li
has no cell outside the spin-1 domain any more (K2: 3 → 0); the window-level
tensor asymmetry moves ×0.984 / ×1.005 / ×1.035 (configs 0 / 1 / 2); HERMES
χ²(b₁) 5.262 → 6.242, χ²(A_zz) 7.308 → 8.191, FAIL → FAIL.

## 1. Every number or label that moved, before → after

Recipes: `scripts/m_kernels.py` (M1, M2, M5 cells, M6, M7 kernel part, the
li6-convolution split), `ab_runs.sh` + `ab_compare.py` + `scripts/m3.py`
(M3), `scripts/runs_extra.sh` + `scripts/cmp_dirs.py` (M4, M9, M10),
`scripts/m5b.py`, the inline M7 run and the M11 timings; all on the shipped
CT18NLO + NNPDFpol1.1 unless said, 300 events seed 11 for the CLI runs.
Declared before any value: the M1 identity at 1e−12 relative, exact equality
at Q² = 2.5 and with the option off.

**The kernel used is the run's.** `shipped_inclusive_kernel(ion,
b1_reference_q2=2.5)` equals the default `Pipeline`'s own kernel in f₁, g₁,
b₁, b₂ at all 30 (x, Q²) points of M1, on d and ⁶Li (bit for bit).

### M1 — b₁/F₁ against Q² (frozen → K1)

Q² = 2.5 / 10 / 100 / 733.6 / 1897 GeV². Under K1 the row is constant; the
largest relative deviation of the unconfined b₁/F₁ from its Q₀² value over
the 60 points is **2.2e−16** (declared 1e−12); at Q² = 2.5 b₁ and b₂ equal the
frozen ones **exactly** (both ions, all six x).

| x | deuteron frozen | deuteron K1 | ⁶Li frozen | ⁶Li K1 |
|---|---|---|---|---|
| 0.1 | 0.00807 / 0.00777 / 0.00721 / 0.00703 / 0.00707 | 0.00807 | 0.002643 / 0.002544 / 0.002361 / 0.002301 / 0.002314 | 0.002643 |
| 0.3 | 1.876e−4 / 2.167e−4 / 2.450e−4 / 2.693e−4 / 2.840e−4 | 1.876e−4 | 6.145e−5 … 9.302e−5 | 6.145e−5 |
| 0.5 | −0.01219 / −0.01605 / −0.02102 / −0.02535 / −0.02776 | −0.01219 | −0.003994 … −0.009091 | −0.003994 |
| 0.7 | 0.03756 / 0.05686 / 0.08722 / 0.1165 / 0.1329 | 0.03756 | 0.0123 / 0.01862 / 0.02857 / 0.03815 / 0.04352 | 0.0123 |
| 0.871 | 0.2548 / 0.2592ᶜ / 0.394ᶜ / 0.473ᶜ / 0.482ᶜ | 0.2548 | 0.08344 / 0.1523 / 0.2897 / 0.4445 / 0.5368 | 0.08344 |
| 0.955 | 0.2592ᶜ at every Q² (unconfined 0.544 … 6.99) | **0.2592ᶜ** (unconfined 0.5442 at every Q²) | 0.1782 / 0.4041 / 0.9854 / 1.776 / 2.291 | 0.1782 |

ᶜ = the J1 note's confinement binds (the value is its upper bound). On the
deuteron at x = 0.955 the table itself (b₁/F₁ = 0.5442 at Q₀²) is above the
cap 3(1 − |g₁|/F₁) = 0.2592, so the confinement still binds there at every
Q² (M5); at x = 0.871 it no longer binds (0.2548 < its bound).

### M2 — the factor F₁(x, Q²)/F₁(x, 2.5) over the default windows

Per x row, [smallest, largest] over the row's accepted cells; the deuteron's
and ⁶Li's are the same to the printed digits (per-nucleon ratio).

| x row | config 0 | config 1 | config 2 |
|---|---|---|---|
| ≈ 0.001 | 0.797 (one Q²) | [0.796, 1.166] | [0.796, 1.748] |
| ≈ 0.01 | [0.839, 1.383] | [0.839, 1.936] | [0.839, 2.283] |
| ≈ 0.1 | [0.978, 1.123] | [0.980, 1.166] | [1.019, 1.170] |
| 0.288 | [0.754, 0.978] | [0.697, 0.935] | [0.677, 0.861] |
| 0.501 | [0.510, 0.791] | [0.438, 0.777] | [0.438, 0.685] |
| 0.724 | [0.319, 0.548] | [0.265, 0.548] | [0.265, 0.504] |
| 0.871 | [0.192, 0.336] | [0.155, 0.336] | [0.155, 0.336] |
| 0.955 | [0.100, 0.154] | [0.078, 0.154] | [0.078, 0.154] |
| whole window | 0.100 (x 0.955) … 1.463 (x 0.0166) | 0.078 … 1.969 (x 0.0079) | 0.078 … 2.367 (x 0.0060) |

So K1 lowers b₁ at high x by up to a factor 13 at the top cell and raises it
at low x and high Q² by up to ×2.4; cells below Q² = 2.5 (the frozen-edge
region of CT18NLO included) get factors 0.80–1.

### M3 — what moved on the runs (A/B, pre `a1f2b4d`, post this tree)

`ab_runs.sh`, 22 CLI runs + 6 examples + 24 inclusive-d (x_max 0.85) + 24
inclusive-⁶Li plan × config runs; exits identical pre and post (only
`inc6_hul` exits 1, A17).

| runs | npz arrays | npz meta | HepMC3 |
|---|---|---|---|
| 17: `coh*` (5), `inc7`, `t6*` (5), `t7*` (3), `tdp*` (3) | identical | `b1_reference_q2` only (absent → its not-read label); knob row `b1_reference_q2` added | identical |
| 6 examples | — | — | identical except `ex_full` (inclusive 6Li through a `Pipeline`): its `GenCrossSection` lines only (496235.884 → 496235.431 pb) |
| 41 of the 48 inclusive + `inc6`, `inc6_rc`, `inc6_vmc`, `incd` | moved (the (x, Q²) draw of a few events per 300: 1–11) | `b1_reference_q2` (2.5), `rank2_input` (+ the K1 clause), `sigma_pb`, `sigma_gen_mb`, `sigma_per_category_pb`; `incd` also `b1_clipped_rate_fraction` | moved (`GenCrossSection`, particle lines; `incd` also its weight line) |
| 10 of the 48 (`d_cos2phi_c1`, `d_transverse-tensor_c1`, `li6_apar_c0/c1`, `li6_cos2phi_c1`, `li6_flip_c1`, `li6_helicity-flip_c0/c1`, `li6_tensor-flip_c1`, `li6_transverse-tensor_c1`) | identical (no event changed cell at 300 events) | as the row above | `GenCrossSection` lines only |

No `lipolgen.*` run attribute changed in any file. The per-category cross
sections (post/pre) and the window-level tensor asymmetry (σ_last − σ_first)/
(σ_last + σ_first) of the plan's categories (`tensor-thirds`, `azz`: azz0 vs
azz+; `tensor-flip`, `flip`: flip0 vs flip+):

| ion (x_max) | config | σ per category, post/pre | tensor asymmetry pre → post |
|---|---|---|---|
| d (0.85) | 0 | azz± 1.0000137, azz0 0.9999726 | 1.293907e−3 → 1.273375e−3 (**×0.9841**) |
| d (0.85) | 1 | 0.9999970 / 1.0000060 | 9.028245e−4 → 9.073303e−4 (**×1.0050**) |
| d (0.85) | 2 | 0.9999861 / 1.0000277 | 6.004660e−4 → 6.212328e−4 (**×1.0346**) |
| d (1.0, `incd`) | 1 | 0.9999970 / 1.0000060 | 9.028669e−4 → 9.073578e−4 (×1.0050) |
| ⁶Li (1.0) | 0 | 1.0000045 / 0.9999910 | 4.240754e−4 → 4.173072e−4 (**×0.9840**) |
| ⁶Li | 1 | 0.9999991 / 1.0000018 | 2.969628e−4 → 2.983333e−4 (**×1.0046**) |
| ⁶Li | 2 | 0.9999955 / 1.0000090 | 1.967244e−4 → 2.035061e−4 (**×1.0345**) |

The `tensor-thirds` total σ is unchanged to 1e−9; the vector plans
(`helicity-flip` = `apar`) move σ by ≤ 9.4e−6 and their asymmetry by ≤ 1e−5
relative (b₁ enters through P_zz of the fill); `transverse-tensor` = `cos2phi`
moves σ by ≤ 6.9e−6. The window-level asymmetry is a rate-weighted average:
low x (factor > 1 at config 2's high Q²) and high x (factor < 1) partly
cancel, hence the per-mille-to-3.5 % shifts despite factors of 0.08–2.4 cell
by cell.

### M4 — the switch, and the deuteron at x_max 1.0

`--b1-ref-q2 0` on post against pre, the same 48 inclusive runs: **48 / 48
npz arrays and 48 / 48 HepMC3 files byte-identical**; the meta differs by
exactly `b1_reference_q2` (absent → 0.0) and the `rank2_input` clause
("--b1-ref-q2 0, the frozen table"). The deuteron at the default x_max 1.0,
all 8 plans × 3 configs: **exit 0 on 24 / 24 pre and 24 / 24 post**.

### M5 — the deuteron's confinement after K1 (all 8 plan names give the same; x_max 1.0)

| config | `b1_clipped_rate_fraction` pre → post | cells pre → post | rows | half |
|---|---|---|---|---|
| 0 | 4.791951543386546e−7 → **1.6772930437472272e−9** | 40 → **15** | 0.871 + 0.955 → **0.955** | upper (joint), 0 lower |
| 1 | 3.871134487837413e−7 → **1.658098263036864e−9** | 59 → **24** | → 0.955 | upper |
| 2 | 3.5361301668316917e−7 → **1.563173218263271e−9** | 59 → **24** | → 0.955 | upper |

The largest unconfined b₁/F₁ is 0.5442 (x = 0.955, at every Q²; banner
"max unconfined b1/F1 = 0.5442 at x = 0.955, Q2 = 206.7", the first cell
reached) against 5.422 / 6.993 / 6.993 before. My own cell-centre recount
reproduces the library's shares to 2e−15 relative. **Is it still needed?** A
deuteron kernel with K1 on and the confinement off (the shipped CT18NLO /
NNPDFpol1.1 objects, MillerB1, the 1e−2 Δ; identical to the shipped one where
nothing binds), as a caller-supplied kernel, x_max 1.0, the CLI fill: **all 24
plan × config runs build** — no refusal by either term. The confinement still
acts (the 15 / 24 cells above) but no default run needs it to build; it stays
(K1 (f)). On the toy pair (`--unpol-sf toy --pol-sf toy`, M9) the deuteron's
share moves 0.0 → 1.26e−8 (ToyF₁ grows with Q² at high x, so K1 raises b₁
there).

### M6 (K2) — ⁶Li after K1

Over each default window (1947 / 3051 / 3498 cells): max b₁/F₁ **0.1782**
(was 1.776 / 2.290 / 2.290), min −0.004739; max(|g₁|/F₁ + b₁/(3F₁)) **0.3421**
(was 0.8746 / 1.0461 / 1.0461); cells outside the spin-1 domain **0** on all
three (was 0 / 3 / 3 at x = 0.955, Q² 1536–1897, 2.372e−11 / 2.712e−11 of the
rate). By K2's own words the item is closed by K1 (no ⁶Li cell of a default
window outside the domain).

### M7 — the price of the reference scale

b₁/F₁ under K1 relative to Q₀² = 2.5 (the same on d and ⁶Li, Q²-independent):

| Q₀² | x = 0.1 | 0.3 | 0.5 | 0.7 | 0.871 | 0.955 |
|---|---|---|---|---|---|---|
| 1.677025 (CT18NLO's floor, `make_lhapdf_unpol_sf(...).q2_min`) | 1.0175 | 0.9480 | 0.9027 | 0.8572 | 0.8003 | 0.7393 |
| 2.5 | 1 | 1 | 1 | 1 | 1 | 1 |
| 4.69 | 0.9807 | 1.0756 | 1.1485 | 1.2313 | 1.3513 | 1.5047 |

Window-level tensor asymmetry (azz0 vs azz+) of inclusive ⁶Li, config 1,
`tensor-thirds`: frozen 2.969628e−4; Q₀² = 1.677025 → 3.287247e−4 (×1.1019 of
the 2.5 value); 2.5 → 2.983333e−4; 4.69 → 2.665315e−4 (×0.8934). The choice
of Q₀² between the grid floor and HERMES's top bin is a ±10 % effect on the
window-level asymmetry and up to −26 % / +50 % on b₁/F₁ at x = 0.955.

### M8 — the HERMES row (`t2_hermes_b1_table2.py`, after line 160 only)

The headline `miller_shipped` is now the run-surface kernel's b₁ at each bin
(`shipped_inclusive_kernel(deuteron(), b1_reference_q2=B1_REFERENCE_Q2)`, the
default `Pipeline`'s own, asserted bit for bit at the six bins in a new test);
the confinement binds at no bin (the confined b₁ equals the unconfined one).
The row's rule, tolerances and other configurations are untouched
(`miller_x1` stays twice the frozen table). Per bin (the three bins below the
grid floor read CT18NLO's frozen-edge F₂; their factor is < 1 through R(Q²)):

| x | ⟨Q²⟩ | F₁(x, ⟨Q²⟩)/F₁(x, 2.5) | b₁ frozen → K1 | pull b₁ | pull A_zz (shipped F₁ᵈ) |
|---|---|---|---|---|---|
| 0.012 | 0.51 | 0.847838 | 5.714659e−2 → 4.845105e−2 | −0.889 → −1.030 | 1.295 → 1.375 |
| 0.032 | 1.06 | 0.905631 | 2.880466e−2 → 2.608640e−2 | −0.837 → −0.924 | 1.108 → 1.169 |
| 0.063 | 1.65 | 0.952691 | 2.092044e−2 → 1.993071e−2 | −1.369 → −1.448 | 1.796 → 1.855 |
| 0.128 | 2.33 | 0.999802 | 8.896023e−3 → 8.894261e−3 | 0.870 → 0.870 | −0.597 → −0.597 |
| 0.248 | 3.11 | 0.980326 | 2.740802e−3 → 2.686879e−3 | −0.043 → −0.058 | 0.131 → 0.143 |
| 0.452 | 4.69 | 0.884239 | −2.066009e−3 → −1.826846e−3 | 1.065 → 1.212 | −0.897 → −1.055 |

| quantity | before (now the frozen sub-row) | after (headline) |
|---|---|---|
| χ²(b₁)/6, bins within | 5.261863562 (4/6) | **6.242226767 (3/6)** |
| χ²(A_zz)/6 on the shipped F₁ᵈ, within | 7.307981761 (3/6) | **8.190725262 (2/6)** |
| χ²(A_zz)/6 on MSTW / R1990, within | 5.549930644 (4/6) | 6.514874394 (3/6) |
| χ²(A_zz)/6 with R1998 | 5.4104 | 6.37872 |
| J4 in-grid sub-row χ²(b₁)/3, χ²(A_zz)/3, within | 1.894165 / 1.177693, 2 / 3 | **2.229766 / 1.490745, 2 / 2** |
| status | FAIL | **FAIL** (the rule's own verdict) |

The frozen sub-row `miller_shipped[frozen table: the reading shipped until K1]
(K1, recorded)` carries the old headline and, nested, J4's old in-grid
numbers. Figure: one extra marker set (grey ×, "Miller x 0.5, frozen table
(until K1): 5.26/6, 4/6 in") and the headline's legend names K1.

### M9 — the toy backend

`ToyF2`'s F₁ depends on Q² (×1.02–1.17 from 2.5 to 10–1897 GeV² at x ≥ 0.5),
so K1 moves the toy runs too: `--unpol-sf toy --pol-sf toy`, inclusive ⁶Li
(config 1, tensor-thirds): arrays moved, tensor asymmetry 2.490430e−4 →
2.511503e−4 (×1.0085); the deuteron: 7.567495e−4 → 7.634065e−4 (×1.0088),
`b1_clipped_rate_fraction` 0.0 → 1.2646e−8. (A test-only consequence: the
toy-named inclusive σ pinned in `test_sf_backend.py` moves 591846.168878409 →
591846.1688784126 pb.)

### M10 — the other ⁶Li b₁ models

`--b1-model li6-convolution --b1-unpol ct18nlo --x-max 0.95`: arrays and
HepMC3 **byte-identical**, meta differs by `b1_reference_q2` = "not read at
b1_model = li6-convolution" only. `--b1-model cdks`: at x_max 0.95 and 0.85
the arrays are identical at 300 events, the category σ move (tensor
asymmetry 8.353213e−6 → 7.060243e−6, ×0.845, at 0.95; ×0.856 at 0.85); at
x_max 1.0 exit 1 pre and post (the rank-2 refusal at x = 0.955). Its band
×2 row (`--b1-band-scale 2`) on the shipped F₁: refused at x_max 0.95 and 0.9
pre (x = 0.871) → **builds at 0.95 and 0.9** post, still refused at 1.0
(x = 0.955).

### M11 — one quick timing (same machine, load ≈ 1, alternating)

`lipolgen-run --channel inclusive --isotope 6Li --events 2000000 --quiet` (no
output), three alternating runs: pre 5.63 / 6.56 / 6.38 s, post 6.19 / 6.36 /
7.19 s (medians 6.38 / 6.36) — no difference visible at this spread. The
example the USAGE §8 recipe times (`generate_inclusive --events 2000000 --out
""`) builds its own class-level kernel and does not read K1: 0.391 → 0.406 M
ev/s, once each (noise). The extra F₁ evaluation is one `NuclearF2::f1a` per
`tables()` call where the option is on.

### Labels, refusals, knob row

- `meta["b1_reference_q2"]`: the number where read, the not-read label
  elsewhere ("not read on channel coherent-6Li", "… on 7Li (no rank-2
  input)", "… at b1_model = li6-convolution", "… on channel tagged-6Li-alpha
  without --inclusive-b1", "not read with a caller-supplied kernel"). **npz
  only**: `docs/HEPMC3_CONVENTION.md` enumerates the `lipolgen.*` keys ("The
  keys are …; One function lists them — `run_provenance`") and may not be
  edited, so the key is not added to `run_provenance` (the
  `b1_clipped_rate_fraction` precedent); T6 is 25 / 25 with the document's
  sha256 unchanged. Said in USAGE §2a where the key is defined.
- Knob row `b1_reference_q2` / `--b1-ref-q2` (70 rows on the default
  `lipolgen-run`, was 69); banner "b1 reference scale (K1 2026-11-17,
  --b1-ref-q2): b1/F1 at Q0^2 = 2.5 GeV2, …"; the `rank2_input` clause.
- Refused: negative / NaN / inf; with CT18NLO 0 < value < 1.677025 (grid
  range [1.677025, 1e10]); any value but 2.5 where not read;
  `struck.b1_reference_q2 ≠ 0`. Accepted: 0, any value in the grid, any
  positive value on a backend that reports no grid (ToyF₂, MSTW, custom).
- **New**: tagged-6Li-alpha with `--inclusive-b1` at the default x_max: pre
  **exit 1** (the struck deuteron's unconfined MillerB1 slot refused at
  x = 0.955, Q² = 167.3, −1.13 from the rank-2 term) → post **exit 0**; at
  x_max 0.85 both build and the tensor asymmetry moves ×1.0047.

### Scope check (`toy_b1(`, `MillerB1`, `CdksB1`, `b1_convolution(` in src/, include/, python/, examples/, validation/)

| site | what it reads now |
|---|---|
| src/core/pipeline.cpp `default_inclusive_kernel`, deuteron branch (`MillerB1`), Miller branch (`Li6B1(MillerB1)`, 6Li), Cdks branch (`Li6B1(CdksB1)`) | run surface: **K1** (the `Pipeline` passes `b1_reference_q2_applied(cfg)`); class level / every factory call: frozen, unchanged |
| src/core/pipeline.cpp `struck_cluster_kernel` (`toy_b1` under `inclusive_b1`) | tagged-6Li-alpha through a `Pipeline`: **K1**; a direct caller: frozen (option default 0) |
| src/core/b1_nuclear.cpp `RawCdksB1` / `cdks_b1_raw_per_nucleon` (inside `Li6ConvolutionB1`) | frozen, unchanged (recorded, §5 item 1) |
| include/lipolgen/sf.hpp `toy_b1`, `b1_convolution`, `MillerB1`, `CdksB1`, `Li6B1`; python/bindings.cpp their classes | the frozen tables themselves (scaled only by the kernel option); unchanged |
| examples/generate_inclusive.cpp (own `Li6B1(MillerB1)` kernel, no `Pipeline`) | frozen, unchanged (class level) |
| validation/b1_li6_table.py, dump_b1_default_li6.py, dump_polligen_reference.py | frozen (they name the class-level table on purpose; reference JSONs untouched) |
| validation/benchmarks/t2_hermes_b1_table2.py | headline **K1** (run-surface kernel); `miller_x1`, `cdks_digitized`: the frozen tables, unchanged |
| t5_delta_moment_bound.py, t5_weighted_unweighted.py, t1_r_world.py (`shipped_inclusive_kernel` / `default_inclusive_kernel` of ⁶Li) | frozen (factory default 0); their report rows did not move (§2) |
| t1_nmc_f2d*.py, t1_hermes_f2d.py, t1_g1d_world.py (`default_inclusive_kernel(deuteron())`) | read F₂ / g₁ only; unchanged |

### Tests re-taken (each with `K1 (2026-11-17)`) or passed the frozen reading explicitly

- Re-taken to K1: `test_bench_spin_deuteron.py` (HERMES headline pins
  6.242 / 8.191, within 3 / 2, MSTW 6.515 / 3, in-grid 2.229766 / 1.490745 /
  (2, 2), `miller_x1 == 2 × b1_frozen_table`, the sub-row list);
  `test_d7_shipped_kernel.py::test_it_is_the_default_pipelines_own_kernel`
  (the Pipeline's kernel is `shipped_inclusive_kernel(li6(),
  b1_reference_q2=2.5)`; the six-point pins stay on the frozen factory);
  `test_j1_deuteron_b1.py` (the slot is MillerB1 × the factor, then confined;
  `b1_clipped_rate_fraction` 1.658098263036864e−9 on 24 cells);
  `test_price_defaults.py` (A_zz −6.599826e−4 → −6.630289e−4, toy −5.534748e−4
  → −5.581584e−4, the four A2 of each, the A_zz ratios mstw 1.0570277 →
  1.0630107, toy 0.8386202 → 0.8418312); `test_sf_backend.py` SIGMA_TOY /
  SIGMA_D7 inclusive (ulp-level: 591846.168878409 → 591846.1688784126,
  496334.15595404804 → 496334.15595404815).
- Frozen passed explicitly (the test pins pre-K1 physics):
  `test_b1_model.py::test_d7_the_band_x2_row_needs_x_max_085_on_the_default`
  (plus a new K1 test: builds at 0.95 / 0.9, refused at 1.0);
  `test_pzz_mode.py::test_the_measured_price_reproduces` (⁶Li rows, sec. A3's
  numbers); `test_rc_low_x_anchor.py` toy builds (the priced A2 table);
  `test_sf_backend.py::test_the_cli_reports_it_as_a_configuration_error_not_a_traceback`
  (`--b1-ref-q2 0`, the refusal it pins).
- New: C++ "K1: the b1 slot read as b1/F1 at a reference scale" (off = bit
  for bit, the identity at 1e−12, equality at Q₀², b₂ set and unset, the
  F₁(x, Q₀²) ≤ 0 rule, the constructor's refusals, scaling before
  confinement, the factory, the run surface, `validate()`); Python
  `test_k1_b1_ratio.py` (23 cases: default meta and label, `--b1-ref-q2 0` =
  the frozen kernel's arrays bit for bit, every refusal, the not-read labels,
  li6-convolution, tagged `--inclusive-b1`, the CLI flag and banner); the
  HERMES K1 test; two `b1_reference_q2` variants in `test_knob_provenance.py`
  (2.5 → 0 on the default; with `--inclusive-b1`).

## 2. What was verified and how

- M4 is the proof that the knob is the only thing that moved: at
  `--b1-ref-q2 0` the 48 inclusive runs are byte-identical to pre in arrays
  and HepMC3. M3: 7Li, coherent, every tagged run at its defaults, the five
  examples without an inclusive `Pipeline`: arrays and HepMC3 identical, the
  meta differing by the new key only.
- The report: `make_report.py --all --plots --compare` against the pre copy
  of `REPORT.json`: **3 changes, all on `t2_hermes_b1_table2`** (the headline
  text, the `miller_shipped` config, the J4 sub-row; the frozen sub-row new);
  no other row moved — so the t5 / t1 rows that use the factories are
  unchanged, as (d) says. Totals 15 pass / 8 fail / 0 blocked (unchanged).
  The 18 other regenerated figures differed in the footer rows only (pixel
  diff: 8 rows at the bottom) and were restored; the HERMES figure stays.
- Line neutrality: `git diff --numstat` is +n/−n on every pre-existing file
  except those that only gained lines at their end (constants.hpp,
  pipeline.hpp, pipeline.cpp, xsec.cpp, lhapdf_sf.cpp, tests/test_pipeline.cpp,
  the end-appended tests) and the HERMES harness (edited after line 160 only;
  its one citation is lines 15–17). Two of my first placements sat in cited
  blocks (the rank-2 dispatch block of `InclusiveKernel::tables` that
  SPIN32_FINITE_GAMMA.md cites, and the `scenario.validate()` use site of
  `PipelineConfig::validate`) and were moved out; the one re-recorded range
  is USAGE §2a's (the range the `--b1-model` wiring row of PHYSICS_CHANNELS
  cites, which K1's clause keeps true).
  The strict and range-audit tallies are identical to `a1f2b4d`'s (measured on
  an export of `a1f2b4d`).
- `li6-convolution` (K1 (c)): confirmed in `b1_nuclear.cpp` that its α–d
  D-wave terms convolve F₁ at the running Q², and found that its embedded-S
  and CG terms do not (§5 item 1); nothing changed there.

## 3. Files

include/lipolgen/constants.hpp, xsec.hpp, sf.hpp, pipeline.hpp;
src/core/xsec.cpp, pipeline.cpp; src/lhapdf/lhapdf_sf.cpp;
python/bindings.cpp, python/lipolgen/__init__.py, cli.py;
tests/test_pipeline.cpp; python/tests/test_k1_b1_ratio.py (new),
test_bench_spin_deuteron.py, test_b1_model.py, test_d7_shipped_kernel.py,
test_j1_deuteron_b1.py, test_knob_provenance.py, test_price_defaults.py,
test_pzz_mode.py, test_rc_low_x_anchor.py, test_sf_backend.py;
validation/benchmarks/t2_hermes_b1_table2.py,
validation/physics_channels_ranges.json (one fingerprint);
docs/USAGE.md (§2a: the two model rows, the K1 paragraph with the knob, the
key and its npz-only status), docs/PHYSICS_CHANNELS.md (the Inclusive row,
the `--b1-model` wiring row); docs/benchmarking/REPORT.json, REPORT.md,
figures/t2_hermes_b1_table2.png; this record; the P1 row of `STATUS.md`.

## 4. Gates (final tree)

Foreground, one long job at a time, 2026-11-19 (machine quiet, load ≈ 1):

| gate | result | wall |
|---|---|---|
| build (`cmake --build build -j2`) | 0 warnings (last full rebuild after the header changes 163 s; the final call up to date) | 31 s last incremental |
| `lipolgen_tests` (full) | **442 / 442 cases, 17 245 521 assertions**, 0 failed (`a1f2b4d`: 441 / 17 244 832; +1 case, K1's) | 687 s |
| pytest, six parts (`pytest_part.sh`) | A 361 / 4 s; B 138 / 1 s; C 188; D 256 / 90 s; E 284 / 61 s; F 347 — **sum 1574 passed, 156 skipped, 0 failed** (`a1f2b4d`: 1535 / 156). Part C first failed 3 doc-link tests because this record then cited two `path:line`s; removed, part C re-run: 188 passed | 209 / 115 / 97 / 417 / 305 / 343 s |
| report `--all --plots --compare` (pre copy) | **15 pass / 8 fail / 0 blocked; 3 changes, all on t2_hermes_b1_table2** (exit 1 by the M2 rule, every change predicted by K1); 18 figures footer-only, restored; the HERMES figure, REPORT.json and REPORT.md stay | 142 s |
| report `--check` | OK | — |
| bench subset (`LIPOLGEN_BENCH=1 … -k bench`), after the report | **373 passed, 1 skipped** (`a1f2b4d`: 372 / 1; +1 the HERMES K1 test) | 493 s |
| A/B (M3, M4) | §1 M3 / M4: the 48 `--b1-ref-q2 0` runs byte-identical to pre | 90 s post set |
| LHAPDF-off build (`build_off.sh`) | configure 0; build 0 warnings; **423 / 423 cases, 17 244 124 assertions** (`a1f2b4d`: 422 / 17 243 435) | 1 + 167 + 467 s |
| citations strict | 1243 / 95 / 7 / **0 broken** / 24 (SPIN32 19 / 115, 0 broken; PYTHIA_BRIDGE 0 broken) — identical to `a1f2b4d` | 3 s |
| citations `--records-only` | 1432 citations in 103 records (+1, this one), **0 broken** | < 1 s |
| citations `--audit-ranges` | 209, **0 REFUSED** (39 by quoted symbol / 103 by exemption, as on `a1f2b4d`); one range re-recorded (USAGE §2a's) | 3 s |
| SPDX | 136 files, all carry the line (the new test and this record carry it too) | — |
| `test_doc_link_gate.py` | **112 passed** | 33 s |
| `git diff --check` | clean | — |
| `git diff --numstat` | +n/−n everywhere except end-appended files (§2) | — |
| `sha256sum docs/HEPMC3_CONVENTION.md` | 2727952ae5a24ec59ca5f6eaf81a159d65517ab2ab0fd3a8c230f36634964eb8 (unchanged) | — |

## 5. Open points

1. **`li6-convolution` is only partly at the running scale (recorded, not
   changed — K1 (c)).** In `Li6ConvolutionB1::terms` (src/core/b1_nuclear.cpp)
   the α–d D-wave terms `dd_d`, `dd_a` convolve F₁ᵈ, F₁ᵅ at the running Q²,
   but the embedded-deuteron S term `s` and the CG D-wave term `cg` convolve
   `b1_d_`, which by default is `cdks_b1_raw_per_nucleon()` — the **frozen
   digitized CDKS Q² = 2.5 GeV² column** (`RawCdksB1`, `kB1CdksQ2p5`) —
   evaluated at the running Q² but Q²-independent. Share of the frozen part
   (s + cg)/total on the default object (`Li6ConvolutionB1()`, ToyF₂ for the
   D-wave terms): (x, Q²) = (0.1, 2.5): 1.68; (0.1, 100): −0.83 (total
   changes sign: s = −4.69e−5 fixed, D-wave 1.91e−5 → 1.04e−4); (0.3, 10):
   1.15; (0.5, 100): 0.72; (0.7, 1000): 0.91. So the dominant term of the
   convolution is a frozen table; reading it as b₁/F₁ at 2.5 would be K1
   applied inside the convolution — the supervisor's to decide.
2. **What P2 must re-take** (stale numbers noticed, not edited here): USAGE
   §7c "69 knobs" (now 70) and its banner example; USAGE §2a / the CLI
   `--b1-band-scale` help and `--x-max` help "the band's x2 row needs about
   0.85" (now builds at 0.95 on cdks); the channel ranking's A_zz price and
   every rank-2 number of ⁶Li / d in USAGE and PHYSICS_CHANNELS (×0.984 /
   ×1.005 / ×1.035 window asymmetry); the J1-note numbers quoted in USAGE line
   162 and in `docs/OPEN_ITEMS_SOLUTIONS.md` §14 (4.79e−7 … , 40 / 59 / 59
   cells, "0.2592 at 0.955, where the unconfined toy's is 5.42") — now
   1.68e−9 / 1.66e−9 / 1.56e−9 on 15 / 24 / 24 cells, unconfined max 0.544;
   USAGE's ⁶Li joint-domain sentence (3 cells → 0, K2); the HERMES numbers
   wherever quoted (5.26 / 7.31 → 6.24 / 8.19); the cli.py comment block on
   the band (D7) numbers; the sampler refusal message's "the mandatory b1
   band's x2 row needs about 0.85". The throughput (J9) needs no re-take
   by M11's reading, but P2's three-run measurement decides.
3. **HepMC3 provenance gap.** A HepMC3 file does not tell the K1 reading
   from the frozen one (`lipolgen.b1_model` = `miller` in both; the new key is
   npz-only because the convention document lists its keys and is frozen).
   If the convention is ever amended, the key belongs beside `b1_model`.
4. Harnesses that call the shipped *factory* (`t5_weighted_unweighted`,
   `t5_delta_moment_bound`, `t1_r_world`) still read the frozen b₁ by
   decision (d); none of their rows moved, but `t5_weighted_unweighted`
   describes its kernel as "the shipped" one — after K1 its b₁ is not the
   run's. A one-word re-label or `b1_reference_q2=2.5` there is P2's / the
   supervisor's call (it would move that row).
5. The toy pair (`--unpol-sf toy --pol-sf toy`) is read by K1 (ToyF₂'s F₁
   evolves), so "naming the toy reproduces the pre-D7 default bit for bit"
   now holds only with `--b1-ref-q2 0`.
6. The default deuteron needs the confinement no more to build (M5), but it
   still binds on 15–24 cells at x = 0.955, where MillerB1's tapered table
   itself has b₁/F₁ = 0.544 > 0.2592 at Q₀² — the table's own top cell, not
   a scale effect.
