<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-26 — phase C5: Boeglin's d(e,e′p)n reduced cross section against the tree's deuteron n(k)

`PLAN.md` phase C5; `../../benchmarking/02_data_nucleon_deuteron.md` §5.4
(and its priority list, item 5: "turns `T25 Hulthén against AV18` from
model-vs-model into model-vs-data out to 550 MeV/c");
`../../benchmarking/08_channel_priority.md` §3.7 (the tagged d–p control,
development order item 1: "d(e,e′p)n momentum distributions (Boeglin), with
the paper's own 10–20 % FSI allowance declared") and §4 ingredient **G** (the
deuteron wave function default). Base `2d89a64`. Owned and written by this
phase: `validation/benchmarks/t2_boeglin_nk.py` (new),
`validation/benchmarks/data/boeglin2011_deuteron_pmiss.csv` (new),
`python/tests/test_bench_boeglin.py` (new), this file, and
`docs/benchmarking/figures/t2_boeglin_nk.png` (written by the harness's
`--plot`). **No physics default, constant, registry label, tolerance of
another row or reference file moved; nothing was committed.** The paper's
PDF stays in the author's corpus (read, never edited); only numbers are
vendored.

## 0. Verdict

| item | verdict |
|---|---|
| the row (`hulthen`, the tagged channel's default) | **PASS (new row).** 22 points of Boeglin *et al.*'s Fig. 2 at θ_nq = 35° (10) and 45° (12), p_m = 0.055–0.538 GeV/c, against s·n(k = p_m), one free scale, the bars ⊕ the paper's own FSI / off-shell allowance A = 20 % ⊕ the digitization: **χ² = 21.98 / 21, p = 0.400, s = 0.481 ± 0.024**. Tolerance p ≥ 0.01, fixed in a hashed design block before any comparison. |
| sub-rows | `av18` **11.54 / 21, p = 0.951**; `cdbonn` **8.89 / 21, p = 0.990** — PASS by the same rule. |
| what decides it | **the allowance.** Without it (bars ⊕ digitization only) every wave fails: 185.9 / 66.5 / 60.3 per 21 (p = 2 × 10⁻²⁸ / 1.3 × 10⁻⁶ / 1.1 × 10⁻⁵). With the paper's own MS FSI factor on the tree's PWIA instead of the allowance, only CD-Bonn passes: 95.6 / 76.7 / **24.4** (p = 1.7 × 10⁻¹¹ / 2.9 × 10⁻⁸ / **0.274**). |
| the default's shape | the worst of the three; d/(s n) reaches 1.25 / 1.41 at 35° and 1.87 / 2.10 at 45° for p_m ≈ 0.50 / 0.54 GeV/c — the Hulthén tail falls too fast, as in the Deeps row. |
| closure of the reading | the paper's CD-Bonn PWIA curve over the tree's CD-Bonn n(p_m) (absolute, no scale): **0.902 (0.876–0.934)** — σ_red = n(p_m) in PWIA is read right to 10 %, the rest the MS off-shell cross section against CC1. |

## 1. The document

W. U. Boeglin *et al.* (Hall A), *Probing the high momentum component of the
deuteron at high Q²*, Phys. Rev. Lett. **107** (2011) 262501,
arXiv:1106.0275v3 — checked against the PDF's first page. The PDF read is the
author's corpus copy `PolarizedLithiumSim/refs/1106.0275.pdf`: 5 pages,
523 354 bytes, sha256 `b9abd413…9c70f91`, pdfTeX output with **no raster
image** on any page. Read with `pdftotext -layout` and `pdftoppm -r 150`,
never with the Read tool.

What it prints: 2H(e,e′p)n at Q² = 3.5 (GeV/c)², beam 5.008 GeV/c; σ_red =
σ_exp/(kσ_cc1) with σ_cc1 the de Forest CC1 off-shell cross section and "the
same form factor as previously (Table I in ref. [30])"; errors: statistical ⊕
the per-bin kinematic systematics, plus "an overall error of 4.5 %" (whether
inside the bars is not said). On the window: "all predict small FSI
contributions (10 - 20%) for 35° ≤ θnq ≤ 45°"; "for 0.2 ≤ pm ≤ 0.5 GeV/c on
average MS predicts between -18% and -5% FSI, JM between -9% and +5% and JVO
between -12% and -16% FSI"; the PWIA calculations "differ at larger pm due to
differences in the wave functions and in details of the off-shell electron
nucleon interaction used"; "within the MS model the calculations using the
CD-Bonn potential are in best agreement with the data". **Fig. 2** (p. 5):
σ_red against p_m at θ_nq = 35°, 45°, 75° (±5°), with the MS model's CD-Bonn
and Paris curves, PWIA and FSI; "The PWIA results are for all angles
identical". **The figure is a vector drawing**: every point, bar and curve is
a drawing operator; the tick labels are glyph outlines (no numeric text on the
page).

## 2. The digitization

**Frames.** Each panel is four 0.64-pt strokes: (a) (114.973, 63.567,
251.737, 258.945) pt, (b) and (c) the same height, 136.765 pt wide each;
p_m = 0 … 0.7 GeV/c, log₁₀ σ_red = −5 … +1. **Calibration check**: the six
major x ticks read back as p_m = 0.10000 … 0.60000 on every panel, and all 55
y ticks (decades and log₁₀ 2 … 9) lie within 2.1 × 10⁻⁵ decade of log₁₀ n; the
labels (0.1 … 0.6; 10¹ at the top, 10⁻⁵ at the bottom) were read from the page
image. The coordinates are float32 (≈ 10⁻⁴ pt), written to 0.001 pt; **0.01 pt
is the declared digitization error** — 5.1 × 10⁻⁵ GeV/c in p_m, 3.07 × 10⁻⁴
decade (0.07 %) in σ_red, negligible against the bars (3–28 %).

**The sibling's tool** (`PolarizedLithiumSim/tools/digitize_figure.py`, read,
never copied): `--inspect` lists the four frames (the panels and the outer
box); with `--frame 114.97 63.57 251.74 258.94 --xrange 0 0.7 --yrange 1e-5
10 --ylog --min-items 5` and one `--curve` per (colour, dash) of the caption
(and the same on panel (b)) it extracts the four curves of panels (a) and (b),
**agreeing with the vendored values to 6.2 × 10⁻⁴** (its frame is rounded to
0.01 pt). It reads no data point or bar, and its label check finds no numeric
word.

**The per-object extraction** is the harness's `digitize()` (dev time,
PyMuPDF 1.28.2): the points are the black 8-arc filled circles in a panel;
each one's bar is the 1.28-pt vertical stroke at its x (one per point,
asserted); the curves are the caption's four (colour, dash) styles, each a
polyline with **exactly one vertex at each data point's p_m** (asserted — the
theory was computed at the bins' averaged kinematics). `--redigitize PDF`
re-reads the PDF and compares the body line for line: **identical, 36 lines**
(also run by the pytest when PyMuPDF and the PDF are present). The bars are
symmetric in σ_red to 1.3 × 10⁻³ of their length.

**The vendored file** `validation/benchmarks/data/boeglin2011_deuteron_pmiss.csv`,
146 lines, sha256 of every byte
`07a6fedecabeadb5124b2dcd3bdce0865e4986151448cad48b13066b9b653334` (written
2026-10-27T00:42:49Z): a provenance block (paper, PDF and sha256, page, figure,
caption, the paper's FSI and error statements, frames, calibration and check,
precision, the commands, the convention — σ_red per deuteron, = n(p_m) with
∫n d³p = 1 in PWIA) and 35 rows — panel, θ_nq, point, p_m, σ_red and its bar's
ends, the four curves at the point, and every PDF coordinate these are read
from:

| panel | θ_nq | points | p_m (GeV/c) | first σ_red (fm³) | last σ_red (fm³) |
|---|---|---|---|---|---|
| (a) | 35° | 10 | 0.17833 … 0.53811 | 0.0148438 ± 0.00068 | 4.898 × 10⁻⁵ ± 0.90 × 10⁻⁵ |
| (b) | 45° | 12 | 0.05499 … 0.53553 | 1.58214 ± 0.047 | 7.679 × 10⁻⁵ ± 2.13 × 10⁻⁵ |
| (c) | 75° | 13 | 0.05656 … 0.53785 | 1.46384 ± 0.031 | 2.876 × 10⁻⁴ ± 0.23 × 10⁻⁴ |

## 3. The design, fixed before the comparison

The verdict's rule, the allowance, the recorded rules (r1)–(r7) and what had
been seen were written into a `==== DESIGN BLOCK ====` of the harness header
and hashed (`437eab47375bd87d96e194c1b5706dd87a52f7b116101776e7adf4920089625d`)
before any comparison was computed; `design_block_sha256()` re-checks the
hash and reads the allowance, the tolerance and the digitization error back
out of the block's sentences on every run (an edited block exits 2, tested).

| UTC 2026-10-27 | event |
|---|---|
| 00:40:35 | the tool's `--inspect` and curve extractions |
| 00:42:49 | `boeglin2011_deuteron_pmiss.csv` written (body from `digitize()`) |
| **00:43:44** | **design block written and hashed** |
| 00:45:08 | first run of the comparisons |
| 00:47–00:50 | post-hoc diagnostic (the TaggedModel cell reading) added outside the block; the figure redrawn through the table's nodes, the 75° points drawn as the block's (r7) says; `--plot` |

**The design.** The 22 points of panels (a) and (b) — the paper's 35°–45°
window up to 0.55 GeV/c; panel (c) is FSI-dominated by the paper's own account
and never compared. Model: n(k = p_m) in fm³ (the paper's own PWIA reduction,
instant form — also how the generator places the spectator) for `hulthen`
(the M-averaged `TaggedModel(deuteron_channel()).n_of_kc`), `av18`
(`source = VmcAV18`) and `cdbonn` ((ψ_s² + ψ_d²)/4π from `cdbonn_wave()`),
times **one free scale common to both panels** ("The PWIA results are for all
angles identical"). **Allowance A = 0.20** — the upper end of the paper's
10–20 %, standing for FSI and the CC1 off-shell prescription together (the
paper quantifies no separate off-shell spread) — as a relative error on each
point, uncorrelated, in quadrature with the bar and the digitization; a
coherent shift is absorbed by s. **PASS iff p(χ², 21) ≥ 0.01**; the default
wave's status is the row's.

**What had been seen before the block was written**, and is said so in it: the
page image; the extraction's printout of all 35 points and the four curves at
each — so the data/PWIA ratios were in view (35°: ≈ 0.65–0.85; 45°: ≈ 0.6–0.9
to 0.46 GeV/c, 1.2–1.3 at the last two points); the `hulthen` and `av18` n(k)
at k = 0.05 … 0.5 GeV and their normalisation integrals; the CD-Bonn
normalisation; the committed Deeps finding. No χ², scale or ratio of the tree
to these data was computed before the block. The allowance and the tolerance
are the task's, fixed before the data were read.

## 4. The numbers (one harness run; window: 35° + 45°, 22 points, p_m 0.055–0.538 GeV/c)

### 4.1 The verdict

| wave | χ² / ndf | p | s (fm³ per n) | max \|pull\| | status |
|---|---|---|---|---|---|
| **`hulthen` (default)** | **21.98 / 21** | **0.400** | 0.481 ± 0.024 | 1.91 | **PASS** |
| `av18` | 11.54 / 21 | 0.951 | 0.567 ± 0.028 | 1.54 | PASS |
| `cdbonn` | 8.89 / 21 | 0.990 | 0.671 ± 0.033 | 1.52 | PASS |

*Corrected in phase F (2026-10-26):* this paragraph first said "The scales
are not comparable across waves (the `TaggedModel` densities are the
generator's piecewise-constant tables, §4.4; CD-Bonn is normalised to 1)".
The cell-step scales above are not directly comparable (the tables
integrate to 1.085 / 1.088, the design block's numbers); read smoothly
(§4.4) all three densities integrate to 1.000, so the smooth-reading scales
0.5015 / 0.5866 / 0.6714 are, and the one free scale absorbs a ×1.34
(Hulthén) / ×1.14 (AV18) excess of n(k) over CD-Bonn in the window — §4.5.
d/(s·n) per point, `hulthen`: 35° 1.06 1.11 0.87 0.89 0.89 0.83 0.87 1.08 1.25
1.41; 45° 1.52 1.57 1.10 0.90 0.88 0.68 0.83 1.01 1.11 1.30 1.87 2.10 —
`cdbonn`: 35° 0.93 1.08 0.91 1.02 0.99 0.88 0.87 1.01 1.07 1.08; 45° 1.19
1.19 1.00 0.88 0.94 0.77 0.90 1.08 1.11 1.20 1.58 1.64.

### 4.2 Recorded (each rule fixed in the block)

- **(r1) absolute closure**: paper CD-Bonn PWIA / tree `cdbonn` n(p_m),
  **mean 0.902, 0.876 … 0.934** over the 22 points; the paper's Paris over
  its CD-Bonn PWIA: 1.161 (0.939 … 1.714) — the wave-function spread the
  paper speaks of, rising with p_m.
- **(r2) theory against theory**, one free scale in ln: the paper's CD-Bonn
  PWIA over s·n is 0.971 … 1.036 for `cdbonn` (the same wave), 0.823 … 1.300
  for `hulthen`, 0.644 … 1.269 for `av18`.
- **(r3) no allowance**: 185.91 / 66.50 / 60.35 per 21 — all fail.
- **(r4) per panel**, with the allowance: 35° alone 4.24 / 7.11 / 1.18 per 9
  (p = 0.895 / 0.626 / 0.999); 45° alone 17.47 / 1.70 / 7.31 per 11
  (p = **0.095** / 0.999 / 0.773).
- **(r5) the MS model's own FSI** in the window: CD-Bonn FSI/PWIA
  **0.746 … 1.109** (mean 0.855) — below 1 − A = 0.8 at six of the 22 points
  (35° at 0.296; 45° at 0.218 … 0.381 GeV/c, down to 0.746), above 1 at the
  two highest 45° points (1.027, 1.109): the paper's figure puts more FSI in
  part of the window than its text's "10 - 20%".
- **(r6) the ratio route**: data / paper CD-Bonn PWIA, 35° mean 0.735
  (0.661 … 0.826), 45° mean 0.837 (0.576 … 1.255); the tree's PWIA times the
  MS FSI factor, A = 0: 95.64 / 76.70 / **24.40** per 21 (p = 1.7 × 10⁻¹¹ /
  2.9 × 10⁻⁸ / **0.274**), s = 0.633 / 0.708 / 0.812.
- **(r7)** the 13 points at 75° are vendored and drawn (grey), never compared.

### 4.3 Where the default misses

Over 0.18–0.42 GeV/c the `hulthen` d/(s·n) is 0.68 … 1.11 (0.68 at 45°,
0.296 GeV/c, where the bar is 20 %); above, it rises: 1.08 / 1.30 at 0.46,
1.25 / 1.87 at 0.50 and 1.41 / 2.10 at 0.54 GeV/c (35° / 45°), against bars of
11–28 % there — it passes because the 20 % allowance and those bars add in
quadrature. Its two lowest 45° points (1.52, 1.57 at 0.055 and 0.0925 GeV/c)
carry the table's cell step (§4.4: +16 % / +8 %). `av18` misses the other way
at 35° (0.81 … 0.71 at p_m ≥ 0.42 GeV/c) and fits 45°; `cdbonn` stays within
0.77 … 1.20 except the two 45° points at 0.50–0.54 GeV/c (1.58, 1.64).

### 4.4 Post hoc: the `TaggedModel` table is piecewise constant

The first figure stepped at low k. That is the declared model:
`TaggedModel.n_of_kc` returns the cell **at or below** k (`tagged.hpp`, the
`TaggedModel` class comment: "the grid IS part of the model"; 280 nodes,
0.0001 … 1.2 GeV, 0.0043 GeV apart), and a cell's value is the density at its
lower node — checked: it equals (ψ̂₀² + ψ̂₂²)/4π from `radial_table` there, to
10⁻¹². At the data's p_m the cell value is **1.005 … 1.160** (`hulthen`) and
1.003 … 1.181 (`av18`) times the same table read smoothly — 1.16 / 1.18 at
p_m = 0.055 GeV/c, 1.08 / 1.09 at 0.0925. A diagnostic, added after the first
run and labelled post hoc in the harness and its output: the verdict's χ² with
the smooth reading is 23.42 / 14.61 / 8.89 (p = 0.32 / 0.84 / 0.99) — the
verdict does not turn on the grid. Beyond this row: the tagged generator's own
spectator density carries the same steps (16–18 % across one cell at
k ≈ 0.05 GeV/c).

### 4.5 Post hoc, phase F (2026-10-26; recorded, not a verdict, outside the hashed block)

Computed by the verification's fix stage with the harness's own `Tree` and
`scaled_chi2` over the same 22 points; the verdict, its allowance and its
tolerance are unchanged.

- **The allowance is the upper end of the paper's estimate.** The paper says
  the models "all predict small FSI contributions (10 - 20%) for
  35° ≤ θnq ≤ 45°" (p. 1); the block's A = 0.20 is its upper end, standing
  for FSI and the CC1 off-shell prescription together. At A = 0.10, the
  paper's lower end: Hulthén 51.57/21 (p = 2.2 × 10⁻⁴, would FAIL), AV18
  27.81 (p = 0.146), CD-Bonn 19.95 (p = 0.525); at A = 0.15: 32.31
  (p = 0.054) / 17.32 / 12.86; at A = 0.20 the recorded 21.98 / 11.54 / 8.89.
  The shipped default passes only in the upper half of the paper's band;
  AV18 and CD-Bonn pass across it.
- **The magnitude the free scale absorbs.** Read smoothly, the three n(k)
  integrate (4πk² n dk over 0–1.19 GeV) to 0.99992 / 0.99983 / 0.99998, so
  the smooth-reading scales 0.5015 / 0.5866 / 0.6714 compare: the free
  scale absorbs a ×1.34 (Hulthén) and ×1.14 (AV18) excess of n(k) over
  CD-Bonn in the window (Hulthén/CD-Bonn 1.24 / 1.50 / 1.40 at k = 0.2 /
  0.3 / 0.4 GeV). With (r1)'s absolute closure (0.876–0.934, mean 0.902)
  and the MS FSI factor (r5, mean 0.855) a right n(k) would need s ≈ 0.77:
  CD-Bonn sits 13 % below that, AV18 24 % and Hulthén 35 %. That would make
  the row an absolute test; it is recorded here, not a verdict — the free
  scale, as declared, discards the magnitude a D2 decision needs.

## 5. What it means for the channels

Ingredient G (`08_channel_priority.md` §4): the deuteron wave function default
moves A_zz^tag by ×2.3–4.8 on the tagged ⁶Li-α channel and sets the d–p
control's spectator spectrum; registry row 9 (D2) keeps the Hulthén pair until
M2. This row is the first model-against-data test of n(k) itself out to
0.55 GeV/c (`02` §5.4's point), and it reads:

- **Under the paper's own allowance the window does not separate the three
  waves** (p = 0.40 / 0.95 / 0.99). The default passes; it is the worst of
  the three, and its misses are where the tagged channels' high-k tail lives.
- **The allowance is what makes room for the Hulthén pair.** On the bars
  alone no pure n(k) shape fits (all three fail at p < 10⁻⁵): FSI and the
  off-shell treatment matter at the level of the bars, as the paper says.
  With the MS model's own FSI factor in place of the allowance, CD-Bonn fits
  (p = 0.27) and Hulthén (p = 2 × 10⁻¹¹) and AV18 (3 × 10⁻⁸) do not — the
  paper's own conclusion, reproduced with the tree's waves, and the same
  ordering as the Deeps row (`t2_deeps_spectator_tail`: CD-Bonn 278.8, AV18
  579.1, Hulthén 627.9).
- For D2 at M2: the Hulthén default is now **η +17σ** (`t2_deuteron_static`),
  **the worst n(k) shape on Deeps and on Boeglin**, and passes here only
  inside the paper's FSI allowance.

## 6. What this does not validate

The tagged spectrum at EIC kinematics (fixed target, Q² = 3.5 GeV², the
neutron unobserved, the target rest frame); the D state separately (n sums S
and D); polarization or the tensor structure (A_zz^tag); the light-cone versus
instant-form question (the instant form is the paper's reduction); FSI itself
(an allowance — `GlauberFsiWeight` is not used); the absolute magnitude of n(k)
(one free scale; the absolute closure is (r1), and the data sit at 0.73 /
0.84 of the paper's own PWIA); p_m above 0.55 GeV/c.

## 7. Priced, not taken (for `DECISIONS_PRICED.md`)

- **Registry row 9 (D2, the deuteron wave default).** Measured here, on the
  22-point 35°–45° window: the Hulthén default χ² 21.98 / 21 (p = 0.40) under
  the paper's 20 % allowance, against AV18 11.54 and CD-Bonn 8.89; with the MS
  FSI factor instead, 95.6 against 76.7 and 24.4 (only CD-Bonn passes). The
  cost of a switch is not measured here (it is the price table's, phase A3).
- **The `TaggedModel` grid** (not a registry row today): the nearest-cell-
  below lookup puts +16–18 % on the density at k ≈ 0.05 GeV/c relative to the
  smoothly read table; interpolating is a code change to a model the header
  says "IS" the grid — listed, not priced.

## 8. Files and commands

Written: `validation/benchmarks/t2_boeglin_nk.py`,
`validation/benchmarks/data/boeglin2011_deuteron_pmiss.csv`,
`python/tests/test_bench_boeglin.py`,
`docs/benchmarking/figures/t2_boeglin_nk.png`, this file. Read only: the
sibling's tool and the PDF.

    cd /home/cpeng/Projects/polli/PoliGen && source env.sh
    python3 validation/benchmarks/t2_boeglin_nk.py                  # the row (0.4 s)
    python3 validation/benchmarks/t2_boeglin_nk.py --plot           # + the figure
    python3 validation/benchmarks/t2_boeglin_nk.py \
        --redigitize ../PolarizedLithiumSim/refs/1106.0275.pdf      # dev time: identical, 36 lines
    python3 -m pytest python/tests/test_bench_boeglin.py -q         # 13 passed
    python3 validation/benchmarks/make_report.py --only t2_boeglin_nk \
        --out-dir <scratch>/report_c5                               # 1 pass

`EXPENSIVE = False` (0.42 s). Not yet in `tiers.py`, the README table or
`BENCHMARK_PLAN.md` §4 — the integrator's (the exact rows are in this phase's
return).

*Phase F (2026-10-26, the verification's fix stage), recorded, nothing
changed:* the vendored CSV's Licence line says "no figure, text or table of
the paper is reproduced", while its own header quotes the Fig. 2 caption and two sentences of p. 4 with attribution —
short attributed quotations, not a document body, but the sentence as written
is false. It is to be reworded ("no figure or table is reproduced; the figure
caption and short sentences are quoted with attribution") at the file's next
revision, which re-pins `DATA_SHA256` in `t2_boeglin_nk.py` and the tests; the digests were
not churned for it now.
