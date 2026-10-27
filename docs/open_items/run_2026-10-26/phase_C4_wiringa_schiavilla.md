<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-10-26 — phase C4: Wiringa–Schiavilla's ⁶Li elastic form factors, a theory anchor for the C0 zero

`PLAN.md` phase C4; `../../benchmarking/04_theory.md` T-20 (its §10 item 1);
`../../benchmarking/08_channel_priority.md` §4 ingredient **F** (radiative
corrections) — the ⁶Li elastic form factor's one consumer is the elastic
radiative tail of the inclusive ⁶Li channels, ranks 1 (cos 2φ) and 4 (A_zz).
Base `2d89a64`. Owned and written by this phase:
`validation/benchmarks/t3_li6_ff_wiringa_schiavilla.py` (new),
`validation/benchmarks/data/ws1998_li6_formfactors.csv` (new),
`python/tests/test_bench_li6_ff_theory.py` (new), this file, and
`docs/benchmarking/figures/t3_li6_ff_wiringa_schiavilla.png` (written once by
the harness's `--plot`). **No physics default, constant, registry label,
tolerance of another row or reference file moved; `data/ff/` was not
created; nothing was committed.** The paper's PDF stays in the author's
corpus (read, never edited); only numbers are vendored.

## 0. Verdict

| item | verdict |
|---|---|
| the row (`ho`, shipped) | **PASS (new row).** The shipped `HoSpin1FF` C0 zero, **q₀ = 3.0998 fm⁻¹**, against the one bound WS98's Fig. 1 gives: its IA+MEC C0 line is drawn non-zero up to its last computed vertex **q_last = 2.9190 fm⁻¹**, so its zero lies above that. Rule (the plan's, applied to the zeros the figure allows): q₀ ≥ q_last − (δ + 0.1 fm⁻¹) = **2.8143 fm⁻¹**, δ = one device unit = 0.0047619 fm⁻¹; distance above the bound **+0.2855 fm⁻¹**. **One-sided: the figure does not locate the zero above q_last, so the upper side is not tested.** |
| sub-row `vmc-ft` | PASS by the same rule (no C0 zero below 5 fm⁻¹, so none where WS98's C0 is drawn non-zero). Recorded beside it: its (fc/Z)² is still 7.0 × 10⁻⁷ at q = 2.92 fm⁻¹, where WS98's IA+MEC C0² has fallen under the figure's 10⁻⁷ floor. |
| against the IA line (recorded) | IA's C0 is drawn to q_last = 3.2000 fm⁻¹: the same rule's bound is 3.0952 fm⁻¹ and the shipped zero meets it by **0.0045 fm⁻¹** only. |
| row 5's data band (recorded) | Under this rule a zero at the UVa FB zero (2.6944 fm⁻¹) would contradict WS98's drawn C0; one at Li *et al.*'s minimum (2.8284 fm⁻¹) would not (by 0.0141 fm⁻¹). The theory anchor and row 5's data band pull the tree's zero in opposite directions (§5). |
| closures of the reading | lines vs symbols **HOLD** (11 points, max \|pull\| 0.54); RFY Table I's q **HOLD** (max \|Δq\| 0.0050 fm⁻¹) and its ratio is **constant**, 0.7116; **μ does NOT hold as declared** (0.8952 / 0.8790 against 0.832 / 0.829: the symbols' q quantization, §4.7) — a post-hoc nominal-q diagnostic gives 0.8414 / 0.8268, inside 2 %; **Q does NOT hold** (\|Q\| read 0.4901 fm² against the printed 0.23 — the figure's C2 is the Rosenbluth F_C2/Z times 3/√2 to 0.5 %). |

## 1. The document

R. B. Wiringa and R. Schiavilla, *Microscopic calculation of ⁶Li elastic and
transition form factors*, Phys. Rev. Lett. **81** (1998) 4317,
arXiv:nucl-th/9807037v2 — checked against the PDF's first page (title,
authors, Argonne / Jefferson Lab / Old Dominion). The PDF read is the
author's corpus copy `PolarizedLithiumSim/refs/nucl-th_9807037.pdf`: 11 pages,
152 769 bytes, sha256 `fc51cb93…91d20e8`, Ghostscript 9.22 output. Read with
`pdftotext -layout` (text) and `pdftoppm -r 150` (the page image), never with
the Read tool.

What it prints: AV18 + Urbana IX VMC wave functions, one-body (IA) and
one- plus two-body (IA+MEC) charge and current operators; Fig. 1 (PDF
page 10) F_L²(q) and F_T²(q) of the ground state, the totals as symbols with
Monte Carlo error bars, the C0 and C2 multipoles as lines (dashed IA, solid
IA+MEC), data from Refs. [8, 28–30]; in the text, μ = 0.829 (IA) / 0.832
(IA+MEC) μ_N, Q = −0.23(9) fm², "the two-body contributions … shift the
minimum to lower values of q", "for q ≥ 3 fm⁻¹, the C2 contribution becomes
dominant, and the shoulder seen in the data is entirely due to this
component". No table. **No number for the C0 zero.**

**The figure is a vector drawing** — `pdfimages` finds no raster; every line,
symbol, error bar and data point is a drawing operator — so the data points
were read too (they are second-hand: WS98's rendering of Refs. [8, 28–30]).

## 2. The digitization

**Frames.** Each panel's box is one `qu` (quad) operator: left (104.40,
214.56, 306.00, 585.12) pt, q = 0 … 4 fm⁻¹, log₁₀ F_L² = −7 … 0; right
(306.00, 214.56, 507.60, 585.12) pt, q = 0 … 4 fm⁻¹, log₁₀ F_T² = −10 … −3.
Every coordinate on the page is a multiple of 0.24 pt (1/300 in): **one such
device unit is the declared digitization error** — 0.0047619 fm⁻¹ in q,
0.0045337 decade in F² (1.05 %; 0.52 % in F).

**Calibration check.** The major tick marks read back as q = 0, 1, 1.9952, 3
(left; the q = 2 tick sits one unit low) and 0, 1, 2, 3, 4 (right), and as
log₁₀ F² within 0.0019 of every decade on both panels. Every "10" label word
sits +0.012 … +0.017 decade and every exponent word +0.133 … +0.138 decade
above its tick, with the exponent equal to that tick's decade.

**The sibling's tool** (`PolarizedLithiumSim/tools/digitize_figure.py`,
read, never copied) was run for its calibration printout and for the two C2
lines:

    python3 tools/digitize_figure.py --pdf refs/nucl-th_9807037.pdf --page 10 \
      --frame 104.40 214.56 306.00 585.12 --xrange 0 4 --yrange 1e-7 1 --ylog \
      --min-items 75 --curve "C2_IAMEC:1,0,0:[] 0"
    (and --min-items 82 --curve "C2_IA:0,0,1:[ 4.8 4.8 ] 0"; the right panel
     with --frame 306.00 214.56 507.60 585.12 --yrange 1e-10 1e-3)

It cannot read this figure alone, for four reasons, all recorded in the
vendored file: `--inspect` finds no frame (it looks for four long strokes;
these frames are one quad each); its per-style collection merges a
calculation's C0 and C2 lines (same colour, same dash — `--min-items`
separates the C2 lines from the C0 lines, nothing separates the C0 lines
from them); it reads no symbol or error bar; its y-label check reads "10" and
the exponent as two words. Its two C2 extractions agree with the vendored
rows to 1.1 × 10⁻⁴ (relative; its float32 coordinates against the 0.24-pt-exact
ones).

**The per-path extraction** is the harness's own `digitize()` (dev time,
PyMuPDF 1.28.2), and `--redigitize PDF` re-reads the PDF and compares the
body line for line — **identical, 572 lines** (also run by the pytest when
PyMuPDF and the PDF are present). Its rules: the lines are the stroked
red-solid / blue-dashed paths of ≥ 60 segments in the left panel, C0 the one
starting at the frame's top-left corner (F² = 1 at q = 0), C2 the one
starting on the floor; each symbol series is a [bars][caps][n symbols] run in
drawing order, the k-th bar on the k-th symbol; the left legend (drawn after
the left frame) is excluded. The lines are sampled every 0.04 fm⁻¹ (8–9
device units); where a line meets the floor the plotting program's clip point
is recorded as `clip`, never as a computed vertex.

**The vendored file** `validation/benchmarks/data/ws1998_li6_formfactors.csv`,
728 lines, sha256 of every byte
`c51d31f9e6e2328dcae7818c4bf1f2a3f759780d6ea5da0d7c0715797dee219b` (written
2026-10-27T00:19:51Z): a provenance block (paper, PDF and sha256, page,
figure, caption, frames, calibration and its check, device unit, the commands,
the publishing convention — per nucleus, F_L² = (fc² + F_C2²)/Z² normalised to
1 at q = 0, F_T² = B(Q²)/Z² — and every column), then 503 record rows and 68
slot rows:

| object | rows | first | last |
|---|---|---|---|
| C0 line IA+MEC / IA | 74 + 1 clip / 81 + 1 clip | q = 0, F² = 1 (both) | vertex q = 2.9190, F² = 1.3395 × 10⁻⁷ / q = 3.2000 |
| C2 line IA+MEC / IA | 95 + 3 clips each | q = 0.1190 | q = 4.0000 |
| F_L² symbols IA+MEC / IA | 20 / 20, bars | q = 0.0476 | q = 4.0000 |
| F_T² symbols IA+MEC / IA | 20 / 20, bars | q = 0.0476 | q = 4.0000, F² = 1.35154 × 10⁻⁸ |
| data F_L² Stanford | 39, bars | q = 0.5333 | q = 3.6429 |
| data F_T² Amsterdam / Saskatoon / Stanford | 9 / 15 / 7, bars | | |
| **slot rows** (IA+MEC) | **68**, q = 0.1190 … 2.9190 fm⁻¹ | F_C0 2.93801, F_C2 −0.00161562, F_M1 4.9768 | F_C0 0.00109798, F_C2 0.00259789, F_M1 −0.0502486 |

**The slot rows are `data/ff/li6_elastic.csv`'s columns** (q, F_C0, F_C2,
F_M1, source) in the tree's normalisation (F_c(0) = Z, F_m the Sachs form,
F_C2 the multipole `TabulatedSpin1FF` turns into F_q = 3F_C2/(2√2 η_A)).
`TabulatedSpin1FF` reads four numbers per line and skips any line whose first
four tokens are not all numbers, so it skips every record row (their
form-factor columns are empty) and reads exactly the 68 slot rows — **tested
on a copy in a scratch `LIPOLGEN_DATA_DIR`** (`data/ff/` stays empty; the
copy's fc, fm and fq at the nodes equal the slot columns to 10⁻⁹, and q
outside [0.1190, 2.9190] fm⁻¹ is refused). The slot rule is in the file and in
`slot_rows()`, and the harness re-derives every slot row from the record rows
on every read: one row per computed IA+MEC C0 vertex with q > 0 that has a
computed C2 vertex at the same x; F_C0 = +Z\|F_C0\|; F_C2 = −Z\|F_C2\| before
the C2 floor gap (the sign of WS98's Q) and +Z\|F_C2\| after; F_M1 from the
IA+MEC F_T² symbols through B = (4/3)η(1+η)F_m², signed + below the F_T zero
and − above, interpolated linearly in ln q (the reader's own rule).

## 3. The design, fixed before the comparison

The verdict's rule, every recorded rule and what had been seen were written
into a `==== DESIGN BLOCK ====` of the harness header and hashed
(`a9f1b3d79dc8a771ef2dbbc0c4397ec0aff44290de90048e165d23b647556a0b`) before any
comparison was computed. `design_block_sha256()` re-checks the hash and reads
16 numerical rules back out of the block's sentences against the constants the
code uses, on every run; an edited block exits 2 (tested).

| UTC 2026-10-27 | event |
|---|---|
| 00:10:52 | the tool's calibration runs and C2 extractions |
| 00:19:51 | `ws1998_li6_formfactors.csv` written (body from `digitize()`) |
| **00:21:38** | **design block written and hashed** |
| 00:25:46 | first run of the comparisons |
| 00:27:51 | post-hoc diagnostics (nominal q) added outside the block; addendum (a)–(c) written |
| (after) | `--plot`, the figure written once; the pytest written |

**The figure does not draw the C0 zero.** Its log axis stops at 10⁻⁷ and the
C0 lines end on that floor (IA+MEC after q = 2.9190, IA after 3.2000 fm⁻¹) and
never re-enter the frame up to 4 fm⁻¹, while the C2 lines dip under it at
2.52 … 2.68 fm⁻¹ and come back (their zero). So "the VMC curve's zero" of the
plan exists on the page only as a lower bound, and the block applies the
plan's tolerance (δ + 0.1 fm⁻¹) to the set of zeros the figure allows: PASS iff
q₀_tree ≥ q_last − (δ + 0.1). An edge with no zero below 5 fm⁻¹ meets it. The
upper side is not tested.

**What had been seen before the block was written**, and is said so in it: the
page image; every drawing operator, the line vertices among them, and by hand
from them q_last = 2.919, the floor crossings, the C2 floor gap, a straight
line through the last two IA+MEC C0 vertices (root ≈ 3.02 fm⁻¹) and a
finite-difference quadratic that turns up before zero; on the tree side the
`ho` edge's (fc/Z)² ≈ 4.7 × 10⁻¹⁰ at 2.93 and ≈ 4.9 × 10⁻⁷ at 2.5 fm⁻¹ besides
the committed 3.0998 fm⁻¹ and 2.5710 fm; on the reference side the C2 line's
low-q coefficient read as \|Q\| ≈ 0.48–0.49 fm² (so the Q closure was expected
not to hold), the q = 0.05 F_T² symbol read as μ ≈ 0.84, and the first Stanford
F_T² point against RFY row I.1 read as a ratio ≈ 0.71. The verdict was
therefore foreseeable; the tolerance is the plan's and the task's, fixed
before all of it. What the worker chose after seeing the figure is how the rule
reads a zero the figure does not locate, and the recorded rules (r1)–(r7).

## 4. The numbers (one harness run; windows as stated)

### 4.1 The verdict and its sub-rows

| edge | q₀_tree (fm⁻¹) | bound vs IA+MEC (q_last 2.9190) | status | bound vs IA (q_last 3.2000), recorded |
|---|---|---|---|---|
| `ho` (shipped) | **3.0998** (3.0997711) | 2.8143; +0.2855 above | **PASS** | 3.0952; +0.0045 above — met |
| `vmc-ft` | none below 5 | — | PASS (no zero) | met (no zero) |

### 4.2 (r1) Where the zero might be — information only

Fits in signed F to the last computed IA+MEC C0 vertices: (e1) line through
the last two, root **3.0252**; (e2) line through the last five, **2.9735**;
(e3) quadratic through the last five and (e4) through the last ten: **no root
below 4 fm⁻¹** (the drawn branch decays nearly geometrically, by 0.67–0.76 per
0.04 fm⁻¹). IA: (e1) 3.4138, (e2) 3.3369. The model spread is the answer: the
figure does not fix the zero above q_last.

### 4.3 (r2) The C0 ceiling above q_last

If the C0 line was computed to 4 fm⁻¹ like the C2 lines (an inference), WS98's
IA+MEC C0² is under 10⁻⁷ on (2.919, 4] fm⁻¹. The `ho` edge's maximum (fc/Z)²
there: 6.31 × 10⁻¹⁰ (under); `vmc-ft`'s: **7.00 × 10⁻⁷ at q = 2.920 fm⁻¹ (above)**.

### 4.4 (r3) The slope radius

ln F² = a₁q² + a₂q⁴ over the computed vertices with 0 < q ≤ 1.0 fm⁻¹:
**IA+MEC 2.5699 ± 0.0091 fm** (digitization), the six variants (q_max 0.6 /
1.0 / 1.4, orders 2 / 3) **2.524 … 2.695 fm**; IA 2.5577 fm. The same fit of the
tree at the same q: `ho` 2.5708, `vmc-ft` 2.5645 fm; by Richardson 2.57100 and
2.57099 fm (t3_li6_rms_devries's 2.5710). Uncompared by rule; the WS98 curve and
the tree agree on the slope to well inside the curve's own model spread.

### 4.5 (r4) The C0 shape, one free scale (ln F², σ = one device unit)

| edge vs line | χ²/ndf | s | WS98/(s·tree) at q = 1.00 / 2.00 / 2.52 / q_last |
|---|---|---|---|
| `ho` vs IA+MEC | 1.894 × 10⁶ / 72 | 4.798 | 0.267 / 1.717 / **10.997** / **43.20** |
| `ho` vs IA | 7.564 × 10⁶ / 79 | 14.15 | 0.095 / 0.756 / 7.746 / 1127 |
| `vmc-ft` vs IA+MEC | 7.34 × 10⁴ / 72 | 1.058 | **1.019 / 1.254 / 1.203** / 0.179 |
| `vmc-ft` vs IA | 9.28 × 10⁴ / 79 | 1.492 | 0.754 / 1.155 / 1.773 / 0.901 |

No tolerance (the digitization error is not a theory error). The shipped HO
shape falls far faster than the VMC charge form factor above 2 fm⁻¹; the
`vmc-ft` edge follows WS98's IA+MEC within 25 % up to 2.5 fm⁻¹ and then stays
above it where WS98's C0 heads for its zero.

### 4.6 (r5) F_M1 against fm; (r6) F_C2 against fq

- **F_M1**: the 20 IA+MEC F_T² symbols against s · (4/3)η(1+η)fm²/Z² (q_z = 1.30
  fm⁻¹, b = 1.85 fm): **χ² = 88.91 / 19, p = 5.2 × 10⁻¹¹, s = 1.024 ± 0.016**;
  s = 1: 91.34 / 20. The pulls: −2.4 … −4.4 across the first peak (q = 0.5 …
  1.25 fm⁻¹, the tree too high) and +2.0 … +3.1 across 2.25 … 3.0 fm⁻¹ (the tree
  too low). WS98's F_T zero lies between the symbols at 1.0 and 1.5 fm⁻¹, the
  minimum symbol at 1.2476 (signed + by the slot rule, so the zero is in
  (1.25, 1.5)); the tree's at q_z = 1.30. The normalisation agrees (s ≈ 1); the
  two-parameter shape does not — as `rc.hpp`'s comment on `LI6_FF_FM_QZ_FM`
  anticipated ("a TWO-parameter shape cannot hit WS98's THREE F_T landmarks").
- **F_C2**: \|F_C2^WS / F_C2^tree\| at the ten computed vertices with q ≤ 0.5
  fm⁻¹: **mean 5.872 (5.723 … 6.064)** — WS98's VMC quadrupole against the
  tree's measured Q = −0.0818 fm², times the 3/√2 of §4.7. WS98's C2 zero lies in
  (2.5190, 2.6762) fm⁻¹; the tree's F_q shares F_point, so its C2 zero is its C0
  zero, 3.0998 fm⁻¹ (asserted). WS98's C2 second lobe peaks at 1.08 × 10⁻⁶
  (q = 3.238 fm⁻¹), where the tree's C2² is 4.2 × 10⁻¹⁴.

### 4.7 (r7) Closures of the reading

| closure | rule (declared) | measured | |
|---|---|---|---|
| μ, IA+MEC / IA | \|μ/μ_printed − 1\| ≤ 0.02, two lowest symbols, linear in q² | 0.8952 (+7.6 %) / 0.8790 (+6.0 %) | **does NOT hold** |
| μ at nominal q (post hoc) | same | **0.8414 (+1.1 %) / 0.8268 (−0.3 %)** | holds (diagnostic) |
| Q, IA+MEC | \|c₀/0.23 − 1\| ≤ 0.05, three lowest C2 vertices, Rosenbluth reading | **\|Q\| = 0.4901 fm²** (+113 %); × √2/3 = **0.2311** | **does NOT hold** |
| lines vs symbols | C0² + C2² within 2σ of the F_L² symbol, 0.25 ≤ q ≤ 2.75 | 11 points, max \|pull\| 0.54 | holds |
| RFY Table I, q | \|q_dig − √q²_printed\| ≤ 0.0095 fm⁻¹ | max 0.0050 | holds |
| RFY Table I, ratio | constant iff max/min ≤ 1.042 | **0.7116** (0.7092 … 0.7144) | constant |

- **μ.** The theory F_T² symbols stand on q = 0.05n (n = 1 … 5) and 0.25n fm⁻¹,
  and the device put the odd multiples of 0.05 half a unit low (0.0476 for
  0.05 — 4.8 % of q, 9.8 % in η). At each symbol's nominal q (all 40 lie within
  one unit of it, checked) the closure gives 0.8414 / 0.8268: the B/Z² reading
  of F_T² holds to about 1 %. Post hoc, labelled so; the declared closure stays
  FAIL. At the nominal q the F_M1 χ² is 77.70 / 19 (the mismatch is not the
  quantization).
- **Q.** The figure's C2 line is, at low q, q²·Q/(2Z) with Q = −0.231 fm² — i.e.
  the Rosenbluth / Donnelly–Sick F_C2/Z (q²Q/(3√2 Z), the convention
  `TabulatedSpin1FF` and `rc.hpp`'s multipole conversion use) **times 3/√2 =
  2.121**, to 0.5 %. The paper prints no formula. Either WS98's plotted C2
  carries that factor or its low-q end is not the printed Q; this phase cannot
  tell which. **The slot rows' F_C2 column inherits it** (it is the figure read
  through the Rosenbluth convention), so a promotion must decide it.
- **RFY.** WS98's seven right-panel Stanford points are RFY Table I's (q to
  0.005 fm⁻¹; their digitized relative errors 0.30 / 0.21 / 0.41 / 0.30 / 0.60 /
  0.50 / 0.93 against RFY's 0.30 / 0.20 / 0.41 / 0.30 / 0.60 / 0.50 / 1.00),
  **scaled by 0.7116** relative to RFY's F_M² read as B/Z² with μ₀ = 0.822. That
  reading is RFY's own: Eq. (5) (p. 860, read from the page image at 300 dpi),
  (dσ/dΩ)_M = (e²/2M_pc²)² [1 + sin²(θ/2)]/(η² sin²(θ/2)) ((J+1)/3J) μ₀² F_M²,
  equals σ_M · B · tan²(θ/2) at θ = π with B = (4/3)η_A(1+η_A)G_M², up to
  1 + η_A ≈ 1 + 2 × 10⁻⁴. And WS98's own theory symbols follow B/Z² (the μ
  diagnostic). So WS98 drew RFY's points 29 % below RFY's published
  normalisation; in RFY's normalisation the data at q = 0.85 … 1.39 fm⁻¹ lie
  ≈ 1/0.7116 = 1.405 × above what WS98 drew. **Row 6's free scale, ŝ = 1.419 ±
  0.206 (RFY data over the tree's F_mag, `phase_B2_rfy.md` §0), is the same
  number** — the tree's F_m and WS98's F_T agree in normalisation (s = 1.024
  above). Recorded for row 6, not resolved here (which normalisation of RFY's
  1966 points later analyses use was not read).

## 5. What it means for the channel

Ingredient F (`08_channel_priority.md` §4): the ⁶Li elastic form factor enters
only the elastic radiative tail of the inclusive ⁶Li channels, where §5 item 3
of that file records the C0 zero's effect on ΔA_zz as ≤ 2.2 × 10⁻⁷ — so nothing
here moves a headline number of ranks 1 or 4. What it settles is the reading
of open item Q1 (the HO refit) and of row 5:

- **The theory anchor bounds the zero from below; row 5's data band sits under
  it.** WS98's full calculation keeps F_C0 non-zero to 2.919 fm⁻¹ (IA: 3.200)
  and attributes the measured shoulder at 2.7–3.6 fm⁻¹ to C2 — whose zero WS98
  puts at 2.52–2.68 fm⁻¹ and whose second lobe (≈ 10⁻⁶) is the shoulder. The
  UVa FB zero (2.6944) is a zero of a C0-only fit to F_L² data, i.e. it absorbs
  the C2 structure; by this rule it contradicts WS98's C0. So row 5's FAIL
  (+0.27 fm⁻¹ above its band) and this row's PASS are not in conflict: they test
  different objects, and a refit of `HoSpin1FF` into row 5's band would fail
  here unless it stayed above 2.814 fm⁻¹.
- **The tree's structural gap is C2, not the C0 zero.** `HoSpin1FF`'s F_q
  shares F_point with F_c, so its C2 vanishes at 3.0998 fm⁻¹ and has no second
  lobe (4 × 10⁻¹⁴ where WS98 has 1.1 × 10⁻⁶); its C0 shape falls ×11 too fast at
  2.5 fm⁻¹ after a free scale (`vmc-ft` within 25 %). Both are what the C0
  shape band (`ho` vs `vmc-ft`) and `fq_scale` already price; this row
  measures them against a calculation for the first time.
- **F_m**: the unfitted (q_z, b) = (1.30 fm⁻¹, 1.85 fm) is right in normalisation
  and in its zero's place (WS98 1.25–1.5) and wrong in shape (χ² 88.9 / 19):
  the first peak too high, the second lobe too low — the region rc.hpp says
  "is where the tail lives".

## 6. What this does not validate

Data (rows 5 and 6 are the data level; WS98 is one VMC calculation, 10 %
underbound, its point radius held near 2.43 fm by construction); the zero's
position above q_last; the form factor where the t-peak integrand reaches
(q ≈ 88–278 fm⁻¹, rc.hpp; WS98 stops at 4 fm⁻¹); the RC tail itself (row 7);
the C2 normalisation (the Q closure); ⁷Li (`HoSpin1FF.for_ion` throws). A PASS
says the shipped zero is not contradicted by this calculation's figure, no
more.

## 7. Priced, not taken (for `DECISIONS_PRICED.md`)

- **Shipping the WS98 table** (`data/ff/li6_elastic.csv`; `tests/test_rc.cpp`'s
  "NO 6Li table is shipped" case would flip). Options: (a) keep `HoSpin1FF`
  (status quo); (b) promote the 68 slot rows by one copy. Measured cost of (b)
  as it stands: the support is 0.1190 ≤ q ≤ 2.9190 fm⁻¹ and the reader refuses
  to extrapolate (tested), so every elastic-tail evaluation beyond 2.919 fm⁻¹
  would throw — a promotion needs a continuation law below 0.119 and above
  2.919 fm⁻¹ (the `vmc-ft` edge's exponential continuation is the in-tree
  precedent); the F_C2 column carries the unresolved 3/√2 (§4.7); F_M1 is
  interpolated between 20 Monte Carlo symbols through its zero. Gain: a C0 shape
  within the theory's own and a C2 with a second lobe.
- **Q1 (the HO refit)**: this row adds a constraint any refit must meet — a C0
  zero ≥ 2.814 fm⁻¹ (IA+MEC; ≥ 3.095 against IA).

## 8. Files and commands

Written: `validation/benchmarks/t3_li6_ff_wiringa_schiavilla.py`,
`validation/benchmarks/data/ws1998_li6_formfactors.csv`,
`python/tests/test_bench_li6_ff_theory.py`,
`docs/benchmarking/figures/t3_li6_ff_wiringa_schiavilla.png`, this file.
Read only: `validation/benchmarks/data/rfy1966_li_magnetic_ff.csv` (row 6's,
sha256 `a3ada6b0…4e593aa8`), the sibling's tool and the PDF.

    cd /home/cpeng/Projects/polli/PoliGen && source env.sh
    python3 validation/benchmarks/t3_li6_ff_wiringa_schiavilla.py            # the row (0.5 s)
    python3 validation/benchmarks/t3_li6_ff_wiringa_schiavilla.py --plot     # + the figure
    python3 validation/benchmarks/t3_li6_ff_wiringa_schiavilla.py \
        --redigitize ../PolarizedLithiumSim/refs/nucl-th_9807037.pdf         # dev time: identical, 572 lines
    python3 -m pytest python/tests/test_bench_li6_ff_theory.py -q            # 18 passed
    python3 validation/benchmarks/make_report.py --only t3_li6_ff_wiringa_schiavilla \
        --out-dir <scratch>/report_c4                                        # 1 pass

`EXPENSIVE = False` (0.40–0.49 s). Not yet in `tiers.py`, the README table or
`BENCHMARK_PLAN.md` §4 — the integrator's (the exact rows are in this phase's
return).
