# Phase A — the G3a/G3b gate mechanics, before anything is rerun

Task A2-research of `docs/open_items/run_2026-09-03/PLAN.md` §Phase A.
Read-only survey: this document is the only file this task wrote.

**CORRECTIONS AFTER THE RERUN (task A2, 2026-09-03, same day).** Three
surgical corrections; the argument of this document stands and its
recommendation was applied unchanged. See
`docs/open_items/run_2026-09-03/phase_A_numbers.md` for the measurements.

1. **Every `design_D_b1_li6.md:NNNN` anchor below is now stale.** A2's
   amendment added 61 lines to that file, most of them inside §5.4, so
   anchors at or after `:1280` have moved down by roughly 10 (the G3a/G3b/G3c
   bullets) to 70 (everything after the amendment block). The quoted text is
   still verbatim; only the line numbers moved.
2. **§0.2, §2.6b and §5.8's "exactly one re-pinned number" is wrong: three
   moved.** `test:515` `got.x_min` did move as predicted (0.237224 →
   0.234833), but the two G3c integrals do **not** survive untouched. §2.6a's
   fix — take the trapezoid on the `x >= 0.01` sub-grid — leaves that
   sub-grid starting at **0.011629**, not at 0.01, because
   `linspace(0.001, 1.59, 300)` has no point at 0.01. That endpoint
   relocation is worth **0.35 %** on `test:512` (2.15370e−4 → 2.14607e−4) and
   **0.71 %** on `test:604` (1.07915e−4 → 1.07152e−4), both outside their
   3e−3 tolerance. Both were re-pinned. The conclusion "the grid change is
   cheap" is unaffected.
3. **§2.4's premise does not transfer to MSTW2008 LO.** "A realistic PDF
   drops the low-x zero below the scan floor" is measured true of CT18NLO
   (zero at 0.0098) and measured **false** of CDKS's own MSTW2008 LO, whose
   low-x zero is at **0.0292** on the old grid — above the old scan floor and
   above the old counting window both. G3a therefore passes on the verdict
   row with or without the recommendation. The recommendation was still
   applied, for the reasons §2.5 gives (which do not depend on the PDF) and
   because the CT18NLO row is now pinned in the doctest, where it does need
   it.

Subjects:
`docs/open_items/run_2026-09-02/design_D_b1_li6.md` §5.4 (the criterion),
`tests/test_b1_nuclear.cpp` (the implementation),
`docs/open_items/run_2026-09-02/phase_D_gate.md` (the measured verdict),
`validation/b1_li6_table.py` (the reproducibility handle).

**Measured versus assumed.** Every number in §2.3, §2.6 and §3.4 below was run
in this tree today (`source env.sh`; `python3` against the installed bindings;
`./build/lipolgen_tests`) and the command is quoted next to it. Everything
else is a quotation with a line anchor. Where a claim is an inference and not a
measurement it says so in the sentence. The full doctest suite was **not**
rerun for this task, so the 363/17203863 baseline is quoted from the brief,
not confirmed here; the 19 `b1_nuclear*` cases were rerun and are green.

---

## 0. The four answers in one page

1. **G3a as written is three constraints wearing one clause: a count, a window
   for the count, and a position tolerance per zero.** The design (`:1285-1289`)
   states the window as `[0.02, 1.0]` and the first-zero tolerance as ±0.08
   about 0.0656 — an interval that reaches down to x = 0. The code implements
   both, faithfully (`test:467`, `test:478`), so the narrower one silently wins
   and the position tolerance is never the operative bound at low x. There is a
   third, *undocumented* bound the design never states at all: the scan grid's
   floor, x = 0.01 (`test:444`, `b1_li6_table.py:133`).
2. **The defensible move is to widen, and it takes three edits, not one.**
   Counting window `[0.02, 1.0]` → `(0, 1.0]`; scan floor 0.01 → **0.001**; and
   G3c's ∫b₁ dx held on the x ≥ 0.01 sub-grid, because `integral_b1` is a
   trapezoid of b₁ = (x·b₁)/x and the extension would move it +23 % on 1/x
   weighting alone (§2.6a). The first two are a pair: **either alone leaves
   G3a aborting** (§2.4, measured). The upper edge 1.0 must **not** move —
   there is a real third zero at x = 1.2204 (ToyF2) / 1.1977 (CT18NLO) that
   "exactly two" exists to exclude. Measured payoff: with all three, the
   realistic-PDF curve that breaks G3a today **passes all three of G3a's
   position clauses**, the low-x zero landing at 0.0098, inside its ±0.08 by
   0.024 (§2.6). Cost: **one** re-pinned number, `test:515` `got.x_min`
   0.237224 → 0.234833 (§2.6b).
3. **`validation/b1_li6_table.py` does NOT write `phase_D_gate.md`.** Its own
   module docstring (`:5-7`) says it does; `main()` (`:263-274`) writes exactly
   one file, `docs/open_items/run_2026-09-02/phase_D_regenerated.md`, and only
   under `--write`. That file does not currently exist in the tree. Measured
   full run: **1 min 40 s wall** (§3.4).
4. **Sixteen files carry the ban or a gate-failure statement, not fourteen.**
   `PLAN.md:116-120`'s list of fourteen omits `docs/DEVELOPMENT_PLAN.md`
   (2 lines) and `python/bindings.cpp` (3 docstrings). Full anchor table in §4.

---

## 1. G3's acceptance criterion, clause by clause, against the code

### 1.1 The design, verbatim

`docs/open_items/run_2026-09-02/design_D_b1_li6.md`, §5.4. First the *stated
gate window*, which is a separate sentence three lines earlier and matters for
everything below:

> **:1277-1281**
> **x range.** CDKS's convolution and Umnikov's warning
> ([hep-ph/9605291](https://arxiv.org/abs/hep-ph/9605291), PLB 391 (1997) 177) put the
> validity floor at **x ≳ 0.1**; above x ≈ 0.8 CDKS themselves say the answer is
> dominated by the high-momentum tail of the wave function. **Gate window: x ∈ [0.10,
> 0.80].** Report x < 0.1 and x > 0.8, gate nothing there.

Then the criterion itself:

> **:1283**
> **Acceptance criterion (three parts, all must hold).**
>
> **:1285-1289**
> * **G3a — shape, hard.** The computed x·b₁ has exactly **two** sign changes in
>   [0.02, 1.0], the first **falling** within **Δx = ±0.08** of 0.0656 and the second
>   **rising** within **Δx = ±0.10** of 0.4572, and its maximum in [0.5, 1.0] sits
>   within **Δx = ±0.10** of x = 0.766. Both reference values come from
>   `b1_landmarks_of_table()`, so a re-digitization moves the target automatically.
>
> **:1290-1299**
> * **G3b — magnitude, soft but recorded.** max |x·b₁| over [0.10, 0.80] agrees with
>   **1.0852 × 10⁻³** — the raw-column maximum, compared against the **raw column**,
>   because the per-nucleon question is *resolved* (§9 Q1) and the gate no longer has to
>   absorb it — **within a factor of 2**, and the ratio is written into
>   `docs/OPEN_ITEMS_SOLUTIONS.md` §10 to three digits. The factor 2 is budgeted, and
>   each piece is reported separately so the residual is attributable: AV18 vs CD-Bonn
>   is ~20 % on P_D and more on the tail; MSTW2008 LO vs ToyF2/CT18 is a factor ~1.5 on
>   F₁ at x ≳ 0.6 (measured); the κ = 1 vs κ choice is a factor 1.5–1.7 (§1.3); and the
>   curve is a *figure digitization*. **A factor 2 is not permission to be a factor 2
>   out — it is the width at which the gate stops being able to distinguish those four.**
>
> **:1300-1303**
> * **G3c — Close–Kumano, reported.** ∫b₁ dx from the new code over the same range,
>   next to `close_kumano_integral(true)` = **+4.592 × 10⁻⁴** and
>   `close_kumano_integral(false)` (Miller). Neither is zero; neither is enforced
>   (`sf.hpp`: "Reported, not enforced"). Record both.

And the consequence clause, which is what makes G3b blocking despite the word
"soft":

> **:1355-1360**
> **Escalation.** If after the checklist the peak ratio is still outside a factor of 2,
> **stop and hand back**: publish the layer-0/1/2 results, the checklist answers and
> the prototype table, leave `Li6ConvolutionB1` in the tree behind the flag with a
> `WARNING: A=2 gate not passed` in its header and in `--help`, and put **no ⁶Li number
> in `OPEN_ITEMS_SOLUTIONS.md`**. Item 10 stays open. That is a legitimate and expected
> outcome for a first-mover calculation and is much better than a tuned number.

One more design line is load-bearing for §2, because it is where the number
0.02 came from:

> **:1305-1312**
> **Honest warning to the implementer — read this before starting.** A prototype written
> for this design (AV18 u,w; CDKS Eqs. 16/17/21 exactly as transcribed at κ = 1; the
> library ToyF2 with R = `r_sigma_lt`; target-mass factor on/off) reproduces the
> **shape** … with sign changes at x ≈ 0.02 and ≈ 0.36 (the independent
> re-run for revision 2; the first revision reported 0.03 and 0.28) rather than 0.066
> and 0.457

### 1.2 The code, clause by clause

All in `tests/test_b1_nuclear.cpp`, `TEST_CASE("b1_nuclear T1 gate layer 3: CDKS
Fig. 4 (G3a hard, G3b/G3c recorded)")`, which opens at **:428** and closes at
**:534**.

| # | design clause (anchor) | code (anchor) | verdict |
|---|---|---|---|
| a | reference landmarks come from `b1_landmarks_of_table()`, "a re-digitization moves the target automatically" (`:1288-1289`) | `const B1Landmarks ref = b1_landmarks_of_table();` **:430**, consumed at **:478**, **:479**, **:485**, **:500** | **agrees** — but see §1.3 D-7 |
| b | "exactly **two** sign changes in [0.02, 1.0]" (`:1285-1286`) | filter `if (got.zeros[i] >= 0.02 && got.zeros[i] <= 1.0)` **:467**, then `REQUIRE(z.size() == 2)` **:475** | **agrees literally**; this is the `REQUIRE` that aborts (§2) |
| c | "the first **falling**" (`:1286`) | `CHECK(sl[0] == -1)` **:476** | agrees |
| d | "within **Δx = ±0.08** of 0.0656" (`:1286-1287`) | `CHECK(std::fabs(z[0] - ref.zeros[0]) < 0.08)` **:478** | agrees; `ref.zeros[0]` = 0.065645 |
| e | "the second **rising**" (`:1287`) | `CHECK(sl[1] == +1)` **:477** | agrees |
| f | "within **Δx = ±0.10** of 0.4572" (`:1287`) | `CHECK(std::fabs(z[1] - ref.zeros[1]) < 0.10)` **:479**, restated **:532** | agrees |
| g | "its maximum in [0.5, 1.0] sits within **Δx = ±0.10** of x = 0.766" (`:1287-1288`) | own scan `if (xs[i] < 0.5 \|\| xs[i] > 1.0) continue;` **:480-484**, then `CHECK(std::fabs(xpk - ref.x_max) < 0.10)` **:485** | agrees |
| h | G3b: "max \|x·b₁\| over [0.10, 0.80]" (`:1290`) | `if (xs[i] < 0.10 \|\| xs[i] > 0.80) continue; win = std::max(win, std::fabs(xb1[i]));` **:496-499** | agrees in intent; the linspace does not land on 0.10 or 0.80 exactly (§1.3 D-8) |
| i | G3b: "agrees with 1.0852e−3 … **within a factor of 2**" (`:1290-1293`) | `CHECK_CLOSE(win, 4.77477e-4, 2e-3); CHECK_CLOSE(ratio, 0.439986, 2e-3); CHECK(ratio < 0.5);` **:504-506** | **DISAGREES** — the code pins the failure, see D-4 |
| j | G3b: "the ratio is written into `OPEN_ITEMS_SOLUTIONS.md` §10 to three digits" (§10, "The A = 2 gate, measured") | `MESSAGE(...)` **:501-503**; the document carries it at `OPEN_ITEMS_SOLUTIONS.md:474`, `:428` | agrees |
| k | G3c: "∫b₁ dx … next to `close_kumano_integral(true)` and `(false)`. Record both." (`:1300-1303`) | `MESSAGE` **:509-511**, `CHECK_CLOSE(got.integral_b1, 2.15370e-4, 3e-3)` **:512**, `CHECK_CLOSE(close_kumano_integral(true), ref.integral_b1, 1e-9)` **:513** | agrees |
| l | — (not in the design) | the computed curve's own landmarks pinned: **:515-520** (`x_min` 0.2372, `xb1_min` −3.30912e−5, `x_max` 0.75508, `xb1_max` 4.77477e−4, `z[0]` 0.022113, `z[1]` 0.377357) | code-only regression pins |
| m | — (not in the design) | the margin clauses `CHECK(z[0] > 0.02)` **:530** and `CHECK(z[0] - 0.02 < 0.005)` **:531**; and `CHECK(std::fabs(z[1] - ref.zeros[1]) > 0.05)` **:533** | code-only; **:531** and **:533** are *upper* bounds on agreement — see D-6 |

The reference column is itself pinned before use, at **:433-442**:
`REQUIRE(ref.zeros.size() == 2)`, `CHECK(ref.zero_slope[0] == -1)`,
`CHECK(ref.zero_slope[1] == +1)`, `CHECK_CLOSE(ref.zeros[0], 0.06564, 1e-3)`,
`CHECK_CLOSE(ref.zeros[1], 0.45718, 1e-3)`, `CHECK_CLOSE(ref.x_min, 0.3324, 1e-3)`,
`CHECK_CLOSE(ref.xb1_min, -1.76814e-4, 1e-4)`, `CHECK_CLOSE(ref.x_max, 0.7657, 1e-3)`,
`CHECK_CLOSE(ref.xb1_max, 1.08521e-3, 1e-4)`,
`CHECK_CLOSE(ref.integral_b1, 4.5920e-4, 1e-4)`.

The scan itself: `const std::vector<double> xs = linspace(0.01, 1.59, 300);`
**:444**, evaluated once per point by `xb1_on_gate_threaded` (**:451**, defined
**:104-120**). The zero finder is `b1_landmarks` (`src/core/b1_nuclear.cpp:760-786`):
it sees a zero **only between two consecutive grid points** (`:767-773`) and
places it by linear interpolation; the extrema are grid points, not
interpolated (`:774-781`); `integral_b1` is a trapezoid of `(x·b₁)/x` over the
whole grid (`:782-784`).

### 1.3 Where the code and the design disagree

**D-1 — the design contradicts itself about where G3 may gate at all.**
`:1280-1281` says *"Gate window: x ∈ [0.10, 0.80]. Report x < 0.1 and x > 0.8,
gate nothing there."* G3a then gates a zero at 0.0656 (below 0.10), counts over
[0.02, 1.0] (both edges outside), and gates a peak position in [0.5, 1.0]
(upper edge outside). Only **G3b** honours the stated window. The code
implements G3a's version. Nothing in the tree reconciles the two sentences.
*This is the tension the task asked to be verified: it is real, it is in the
design, and it is wider than the ±0.08 issue alone.*

**D-2 — the counting window is narrower than the position tolerance it
brackets.** `ref.zeros[0]` = 0.065645, so `|z0 − ref.zeros[0]| < 0.08` admits
z0 ∈ (−0.014355, 0.145645), i.e. **any non-negative zero below 0.1456**. The
counting window `[0.02, 1.0]` excludes z0 ∈ [0, 0.02). Both are enforced
(**:467**, **:478**), so a zero at 0.015 passes the stated tolerance and is
nevertheless dropped from the count, taking `z.size()` to 1 and aborting the
case at **:475** before **:478** is ever evaluated. The tolerance is therefore
dead code below x = 0.02.

**D-3 — a third bound exists that the design never states: the scan floor.**
The design fixes no x grid for the computed curve. The doctest chose
`linspace(0.01, 1.59, 300)` (**:444**, and again at **:582** in the checklist
case); `validation/b1_li6_table.py` chose `x_grid(0.01, 1.0, n_x)` (**:133**,
`:136`). A zero below 0.01 is invisible to `b1_landmarks` no matter what the
counting window says. For a realistic PDF that is exactly what happens
(measured, §2.3), so today the *operative* low-x bound is 0.01 and not 0.02.

**D-4 — G3b's code asserts the negation of G3b's criterion.** Design `:1290-1293`
accepts ratio ∈ [0.5, 2]. The code at **:504-506** pins
`CHECK_CLOSE(win, 4.77477e-4, 2e-3)`, `CHECK_CLOSE(ratio, 0.439986, 2e-3)` and
`CHECK(ratio < 0.5)`. That is deliberate and is authorised by the design's own
T1 row (`:1380`: *"G3b/G3c `MESSAGE`d, and `CHECK`ed against whatever the
implementer records, so a later regression is caught"*), and the comment at
**:487-494** says so. The consequence for A2/A7 is mechanical and must be
planned for: **these three lines fail the moment the ratio moves**, in either
direction. The doctest is not the gate; `phase_D_gate.md`'s verdict table is.

**D-5 — "soft" does not mean "not blocking".** `:1290` calls G3b *"soft but
recorded"*, `:1283` says all three parts *"must hold"*, and `:1355-1360` makes
failure of exactly this clause the trigger for Escalation. "Soft" describes the
doctest's treatment, not the gate's.

**D-6 — two code-only clauses assert that agreement is not too good.**
**:531** `CHECK(z[0] - 0.02 < 0.005)` and **:533**
`CHECK(std::fabs(z[1] - ref.zeros[1]) > 0.05)` pin the *margins* as they stand
today. They are honest documentation of a thin pass, but they are upper bounds
on agreement: if MSTW moves the second zero from 0.3774 towards the digitized
0.4572 by more than 0.030 it will trip **:533** by getting *better*. Same class
of edit as D-4.

**D-7 — "a re-digitization moves the target automatically" is only half true.**
The tolerances do float with `b1_landmarks_of_table()` (**:478**, **:479**,
**:485**), but **:436-442** pins that function's output to typed constants at
1e-3/1e-4, so a re-digitization fails those `CHECK_CLOSE`s first. That is a
deliberate regression guard, not a bug; it does mean the automatic-retarget
sentence overstates what happens.

**D-8 — G3b's window edges are not on the doctest grid.**
`linspace(0.01, 1.59, 300)` has step 0.00528428; 0.10 and 0.80 fall at indices
17.03 and 149.5, so the doctest's effective G3b window is
[0.10527, 0.79792] rather than [0.10, 0.80]. `b1_li6_table.py` uses
`x_grid(0.10, 0.80, n_x)` (`:136`), which lands on both edges. This is one of
the two grid differences behind the 4th-digit disagreement that
`phase_D_gate.md:14-19` records; the peak is at 0.755, far from either edge, so
today it changes nothing.

**D-9 — G3a is enforced in exactly one place.** `validation/b1_li6_table.py`
applies **no** counting window and asserts **no** count: `gate_report` returns
`list(lm.zeros)` unfiltered (`:140`) and `emit` prints them (`:234-237`). The
script reports; the doctest gates. When the doctest's window changes, nothing
in the script has to change for consistency — but the script's *grid floor*
does (§2.5).

---

## 2. The counting window — what to do, and why that number

### 2.1 The question, stated precisely

`REQUIRE(z.size() == 2)` at **:475** is a hard abort inside a counting window
whose floor (0.02) is tighter than the position tolerance the same clause
grants (down to 0), and tighter still than the *scan* floor nobody wrote down
(0.01). With a realistic PDF the low-x zero falls below both. So the question is
not "does the physics agree" but "which of three inconsistent bounds is the
gate".

### 2.2 Where the 0.02 came from

The design's prototype table (`:1305-1312`) reports the prototype's own first
sign change at **x ≈ 0.02**, and notes in the same sentence that revision 1 of
the same prototype reported **0.03**. The counting floor is that number. It is
therefore not a validity floor, not a kinematic bound and not a digitization
limit — it is one run's measured answer, known from the design's own text to
have moved by 50 % between two runs of the same code. **A gate floor set at a
quantity the design itself records as unstable at the 50 % level is not a
bound.** By contrast every other number in G3a (0.0656, 0.4572, 0.766) is
computed from the digitized column at run time.

### 2.3 What the curve actually does below x = 0.02 — measured

`DeuteronConvolutionB1` at its default options (`finite_q_delta = true`,
`r1998`, target mass on), Q² = 2.5, x·b₁ on a log-ish probe grid.
Run: `python3 <scratch>/lowx_zero.py` (both δ-functions, 8.0 s wall for the κ
pair, 0.9 s for the κ = 1 pair).

| x | ToyF2, κ = 1 | ToyF2, κ (**default**) | CT18NLO, κ = 1 | CT18NLO, κ (**default**) |
|---|---|---|---|---|
| 0.0001 | +2.674587e−5 | +2.674587e−5 | +8.496117e−6 | +8.496118e−6 |
| 0.0002 | +2.332149e−5 | +2.332149e−5 | +1.005575e−5 | +1.005575e−5 |
| 0.0005 | +1.924132e−5 | +1.924133e−5 | +1.234191e−5 | +1.234192e−5 |
| 0.0010 | +1.636533e−5 | +1.636533e−5 | +1.328640e−5 | +1.328640e−5 |
| 0.0020 | +1.349026e−5 | +1.349030e−5 | +1.232434e−5 | +1.232438e−5 |
| 0.0030 | +1.183292e−5 | +1.183289e−5 | +1.037634e−5 | +1.037634e−5 |
| 0.0050 | +1.254794e−5 | +1.254802e−5 | +8.678571e−6 | +8.678564e−6 |
| 0.0070 | +6.969007e−6 | +6.968637e−6 | +2.556756e−6 | +2.556172e−6 |
| **0.0100** | +5.029761e−6 | +5.029002e−6 | **−6.425156e−7** | **−6.439858e−7** |
| 0.0150 | +2.559705e−6 | +2.557214e−6 | −3.391196e−6 | −3.395498e−6 |
| 0.0200 | +6.634818e−7 | +6.580715e−7 | −4.587856e−6 | −4.596432e−6 |
| 0.0300 | −2.212618e−6 | −2.228564e−6 | −6.480667e−6 | −6.504327e−6 |
| 0.1000 | −1.992742e−5 | −2.033946e−5 | −3.914759e−5 | −4.001790e−5 |
| **falling zero** | **0.022307** | **0.022280** | **0.009398** | **0.009396** |

**The zero is not gone. It exists, it is single, and it sits at x = 0.009396**
at the gate's default δ-function with CT18NLO — 0.0006 below the scan floor of
0.01 and 0.011 below the counting floor of 0.02. Its distance from the
reference is |0.009396 − 0.065645| = **0.056249**, comfortably inside G3a's own
±0.08 (it clears by 0.024, i.e. by more than the second zero clears its ±0.10).
Both curves are positive and monotone-ish over the whole decade [1e−4, 5e−3];
there is no second crossing hiding below the probe.

Also measured, on the doctest's own grid: the ToyF2 default curve has **three**
zeros, `0.022113(falling), 0.377357(rising), 1.220409(falling)`. The third is
what the counting window's *upper* edge of 1.0 exists to exclude.

### 2.4 The window and the floor are a pair; either alone is a no-op

- **Widen the counting window only** (`test:467` floor 0.02 → 0), scan floor
  left at 0.01: the CT18/MSTW zero at 0.0094 is *not on the grid*, so
  `got.zeros` contains one entry ≤ 1.0 and `REQUIRE(z.size() == 2)` still
  aborts. Confirmed by today's script run, whose grid floor is 0.01 and which
  reports `CT18NLO, kappa | 0.43877` — one zero, no counting window applied at
  all (§1.3 D-9).
- **Lower the scan floor only** (0.01 → 0.001), counting window left at
  `[0.02, 1.0]`: the zero at 0.0094 is found and then discarded by the filter
  at **:467** because 0.0094 < 0.02. `z.size()` is 1 and the `REQUIRE` aborts.

So edits (i) and (ii) of §2.5 have to land together. Saying "widen the
counting window" without also lowering the grid floor is the failure mode to
avoid; edit (iii) is a separate correction that neither of them enables, but
that both of them make necessary.

### 2.5 The recommendation, and why these numbers

**(i) Counting window: `[0.02, 1.0]` → `(0, 1.0]`.** Concretely, at
`test:467`, drop the lower comparison and keep `got.zeros[i] <= 1.0`.

Argument. The clause's physics content is the *position* tolerance: "the
computed first zero lands near the digitized first zero". The count and its
window are bookkeeping that makes the position clause well-posed — you must
know which zero is "the first" and you must be able to say that no spurious
extra crossing appeared. Bookkeeping must therefore be **at least as wide** as
the physics it brackets, or it becomes a second, tighter, unstated position
cut, which is exactly what D-2 describes. The widest the design grants is
`ref.zeros[0] − 0.08 = −0.014355`, clipped at the physical floor x > 0.
That gives `(0, 1.0]` and **introduces no new number**: it is the tolerance the
design already wrote down, nothing more.

The **upper** edge stays at 1.0, unchanged. The second zero's tolerance reaches
only `ref.zeros[1] + 0.10 = 0.5572`, so the interval [0.56, 1.0] is doing no
tolerance work — but it is doing real gate work: the measured third zero at
x = 1.2204 shows that the curve does cross again, and "exactly two" is the only
clause that would catch such a crossing migrating below 1.0. Widening the
ceiling would delete a live constraint. Narrowing it to 0.5572 would too.

**(ii) Scan floor: 0.01 → 0.001**, at `test:444`, `test:582` and
`validation/b1_li6_table.py:133` (anchor read 2026-09-15 "floor as part of G3a") (and `:130` if the G3b sub-grid is to stay
inside the union).

Why 0.001 and not another number, in order of weight:

1. **It is a decade below the measured crossing.** CT18NLO at the default
   δ-function crosses at 0.009396 (§2.3); ToyF2 at 0.022280. A floor of 0.001
   is 9.4× below the realistic-PDF crossing and 22× below the toy one, so the
   crossing is bracketed by roughly two decades of grid on the left. MSTW2008
   LO is a real global fit of the same family as CT18; a floor an order of
   magnitude below CT18's answer is the margin that makes the clause robust to
   *which* real PDF is used, which is the whole complaint the docs raise about
   the low-x zero.
2. **The curve is measured to be positive and featureless over
   [1e−4, 5e−3]** for both PDFs (§2.3), so nothing is hidden below 0.001 that a
   still-lower floor would find. This is a measurement, not an assumption.
3. **0.01 is disqualified**: it is the value that fails today, and it is also
   the digitized table's own floor (`design:1259`: *"x ∈ [0.0100, 1.590], 300
   points"*), so keeping it looks principled and is not — the *reference*
   curve's first zero is at 0.0656, an order of magnitude above its own floor,
   and never comes near it. There is no like-with-like argument that requires
   the computed curve to stop where the digitized one stops; the comparison at
   **:478** is against a *number*, 0.065645, not against the table evaluated at
   z0.
4. **1e−4 or lower is disqualified on cost and on honesty.** The design puts the
   convolution's validity floor at x ≳ 0.1 (`:1277-1281`) and
   `phase_D_gate.md:361-363` records the low-x tail as the one thing still
   unexplained. Buying three decades of grid in a region the design says has no
   validity, at 0.24 s per finite-|q| point, is spending the gate's runtime on
   a number nobody will quote.
5. **It keeps the grid a plain `linspace`.** `linspace(0.001, 1.59, 300)` has
   step 0.0053144 against today's 0.0052843 — a 0.6 % change. No new grid
   *shape* (log-spacing, splicing) is introduced, so the one thing that changes
   about every pinned landmark is a small, uniform, explainable shift.

**(iii) G3c's ∫b₁ dx must stay on the x ≥ 0.01 sub-grid.** This is not
optional and it is the one place the grid change would otherwise do real
damage: `integral_b1` is a trapezoid of b₁ = (x·b₁)/x, so extending the grid to
0.001 moves it by **+23 %** for reasons that are pure 1/x weighting. §2.6a gives
the measurement and the design line (`:1300-1301`) that forbids it.

**What this costs — measured, and it is small.** Every landmark the doctest
pins is grid-dependent; `phase_D_gate.md:14-19` makes exactly that point about
the doctest and the script differing in the 4th digit. §2.6b runs all nineteen
pinned quantities on both grids. The result: with (iii) in place, **exactly one
pin needs a new number** — `test:515` `got.x_min`, 0.237224 → **0.234833**,
because the argmin hops one grid point on a flat stretch (its *value*
`xb1_min` moves by 3e−5 relative). Everything else stays inside its existing
tolerance, two of them (`test:517`, `test:603`) with a quarter to a third of the
tolerance to spare. Without (iii), `test:512` and `test:604` break as well.

**The alternative that preserves every pin, and why it is now third choice.**
Keep `linspace(0.01, 1.59, 300)` for everything, and locate the first zero on a
separate small low-x bracket (e.g. `logspace(-3, log10(0.02), 25)`) spliced in
front of it. The union then covers [0.001, 1.59] with an overlap and every
existing pin survives untouched. Before §2.6b this looked like the cheaper
route; it is not — it is two grids to explain and it makes `integral_b1`
ambiguous (over which grid?), to save re-pinning a single number. **Do not**
choose the remaining option — a bare `logspace` for the whole scan — because it
re-pins everything *and* changes the trapezoid weighting of `integral_b1`,
which G3c reports.

### 2.6 The proposed grid, measured end to end

Reproduction of `tests/test_b1_nuclear.cpp`'s layer-3 scan in Python against the
same `DeuteronConvolutionB1` defaults, on the current grid and on the proposed
one. The current-grid ToyF2 row reproduces the doctest's pins exactly
(0.022113 / 0.377357 / 0.755084 / 4.774769e−4 / 2.153703e−4 / 0.439986), which
is what licenses reading the other three rows as the answer the doctest would
give.

Run: `python3 <scratch>/grid_effect.py {toy_cur,toy_new,ct_cur,ct_new}`,
~72 s per finite-q row.

| | ToyF2, grid [0.01, 1.59] (today) | ToyF2, grid [0.001, 1.59] | CT18NLO, [0.01, 1.59] | CT18NLO, [0.001, 1.59] |
|---|---|---|---|---|
| all zeros found | 0.022113(−), 0.377357(+), **1.220409(−)** | 0.022040(−), 0.377374(+), **1.220437(−)** | 0.438771(+), **1.197727(−)** | **0.009795(−)**, 0.438761(+), **1.197722(−)** |
| count in `[0.02, 1.0]` (today's clause) | 2 → **passes** | 2 → passes | **1 → ABORTS** | **1 → ABORTS** |
| count in `(0, 1.0]` (proposed clause) | 2 → **passes** | 2 → **passes** | 1 → ABORTS | **2 → PASSES** |
| \|z0 − 0.065645\|, tol 0.08 | 0.043531 ✓ | 0.043605 ✓ | — | **0.055850 ✓** |
| \|z1 − 0.457177\|, tol 0.10 | 0.079820 ✓ | 0.079803 ✓ | — | **0.018416 ✓** |
| peak in [0.5, 1.0], tol 0.10 on position | 4.774769e−4 @ 0.755084, Δ 0.010579 ✓ | 4.774939e−4 @ 0.755642, Δ 0.010021 ✓ | 7.807455e−4 @ 0.733946, Δ 0.031717 ✓ | **7.807344e−4 @ 0.734385, Δ 0.031278 ✓** |
| G3b `win` / ratio | 4.774769e−4 / **0.439986** | 4.774939e−4 / 0.440001 | 7.807455e−4 / **0.719442** | 7.807344e−4 / 0.719432 |
| `integral_b1` | +2.153703e−4 | **+2.653156e−4** | +2.094970e−4 | **+2.474269e−4** |

**Three readings.**

1. **The proposal does what it is supposed to do.** With the widened counting
   window *and* the lowered floor, CT18NLO — the stand-in for the realistic PDF
   that breaks G3a today — passes all three of G3a's position clauses, and the
   first zero clears its ±0.08 by 0.024. With either edit alone it still aborts
   (columns 3 and 4, row "count in `[0.02, 1.0]`"), which is the claim of §2.4
   measured rather than argued.
2. **The third zero is real on every configuration**, at 1.2204 (ToyF2) and
   1.1977 (CT18NLO). The counting window's ceiling of 1.0 is load-bearing and
   must not move.
3. **`integral_b1` is the one quantity the grid change really moves**: +23 %
   for ToyF2, +18 % for CT18NLO. That is not noise — see §2.6a.

#### 2.6a `integral_b1` must NOT be taken over the widened grid

`b1_landmarks` computes `integral_b1` as a trapezoid of **b₁ = (x·b₁)/x** over
the whole grid it is handed (`src/core/b1_nuclear.cpp:782-784` (as of b1071b1)). Adding
[0.001, 0.01] adds a region where x·b₁ is ~1.6e−5 and 1/x is ~10³, so b₁ itself
is ~1.6e−2 there — two orders of magnitude above anything in the gate window.
The integral is dominated by the extension.

The design already forbids this. G3c is *"∫b₁ dx from the new code **over the
same range**, next to `close_kumano_integral(true)` = +4.592 × 10⁻⁴"*
(`:1300-1301`), and that reference is computed over the digitized table's own
[0.0100, 1.590] (`b1_landmarks_of_table`, `src/core/b1_nuclear.cpp:802-808`).
Letting the computed integral run from 0.001 while the reference starts at 0.01
breaks the like-for-like G3c comparison silently — the number would move 23 %
for a reason that has nothing to do with physics.

**So: lower the scan floor for the zero finding, and keep G3c's trapezoid on
the x ≥ 0.01 sub-grid.** Mechanically, call `b1_landmarks` a second time on the
`x >= 0.01` slice (reusing the already-computed `xb1` values — the expensive
part is the evaluation, not the landmark pass) and take `integral_b1` from
that. Then `test:512`'s pin 2.15370e−4 and `test:604`'s 1.07915e−4 survive
untouched, and G3c keeps comparing the same range to the same range.

#### 2.6b Exactly which pinned assertions move

Measured, ToyF2 at the default δ-function for the layer-3 case (`test:428-534`)
and at κ = 1 for the checklist case (`test:536-606`; run
`python3 <scratch>/checklist_grid.py`, 1.7 s). "current" reproduces every pin
in the tree, which is the check that the harness is faithful.

| pin | quantity | current | on [0.001, 1.59] | rel. Δ | tol | verdict |
|---|---|---|---|---|---|---|
| `test:504` | `win` (G3b) | 4.77477e−4 | 4.77494e−4 | 3.6e−5 | 2e−3 | holds |
| `test:505` | `ratio` (G3b) | 0.439986 | 0.440001 | 3.4e−5 | 2e−3 | holds |
| `test:512` | `got.integral_b1` (G3c) | 2.15370e−4 | 2.65316e−4 | **2.3e−1** | 3e−3 | **BREAKS** — but see §2.6a: it must not move |
| `test:515` | `got.x_min` | 0.237224 | 0.234833 | **1.0e−2** | 1e−3 | **BREAKS** — a pure grid-point relocation, one step |
| `test:516` | `got.xb1_min` | −3.30912e−5 | −3.30902e−5 | 3.2e−5 | 3e−3 | holds |
| `test:517` | `got.x_max` | 0.755084 | 0.755642 | 7.4e−4 | 1e−3 | holds (by 26 %) |
| `test:518` | `got.xb1_max` | 4.77477e−4 | 4.77494e−4 | 3.6e−5 | 2e−3 | holds |
| `test:519` | `z[0]` | 0.022113 | 0.022040 | 3.3e−3 | 5e−3 | holds |
| `test:520` | `z[1]` | 0.377357 | 0.377374 | 4.5e−5 | 2e−3 | holds |
| `test:530` | `z[0] > 0.02` | 0.022113 | 0.022040 | — | — | holds |
| `test:531` | `z[0] − 0.02 < 0.005` | 0.002113 | 0.002040 | — | — | holds |
| `test:593` | `lr.zeros[1]` (κ=1, toy R) | 0.365215 | 0.365202 | 3.6e−5 | 3e−3 | holds |
| `test:594` | `l1.zeros[1]` (κ=1, r1998) | 0.392321 | 0.392320 | 2.6e−6 | 3e−3 | holds |
| `test:595` | `lr.xb1_max` | 2.99437e−4 | 2.99431e−4 | 2.0e−5 | 3e−3 | holds |
| `test:596` | `l1.xb1_max` | 2.94842e−4 | 2.94839e−4 | 1.0e−5 | 3e−3 | holds |
| `test:602` | `l1.zeros[0]` | 0.022137 | 0.022060 | 3.5e−3 | 5e−3 | holds |
| `test:603` | `l1.x_max` | 0.739231 | 0.739699 | 6.3e−4 | 1e−3 | holds (by 37 %) |
| `test:604` | `l1.integral_b1` | 1.07915e−4 | 1.57861e−4 | **4.6e−1** | 3e−3 | **BREAKS** — same cause as `test:512` |
| `test:605` | `l1` G3b ratio | 0.271692 | 0.271689 | 1.1e−5 | 2e−3 | holds |

**So the grid change is cheap.** Of the nineteen pinned quantities across the
two test cases, **two** break for the ∫b₁-range reason that §2.6a says must
be fixed rather than re-pinned, and **one** (`test:515`, `got.x_min`) breaks
because the argmin hops one grid point on a flat part of the curve — its
*value* `xb1_min` moves by 3.2e−5 relative. Re-pin `got.x_min` to **0.234833**;
fix the integral by range rather than by re-pinning. Two pins, `test:517` and
`test:603`, hold with 26 % and 37 % of their tolerance left, so they are worth
a glance in the actual rerun.

**One caveat on the resolution of the low-x zero.** The uniform 300-point grid
brackets CT18NLO's zero between 0.006314 and 0.011629 and interpolates it
linearly to 0.009795; the dedicated log probe of §2.3 puts it at 0.009396. So
the *position* of the low-x zero is only resolved to ≈ 4 % on this grid. That
is irrelevant to the clause — the tolerance is ±0.08 and the discrepancy is
0.0004, five thousandths of it — but it should not be quoted to five digits in
any document.

### 2.7 Why "redefine the clause" is the wrong move here

For completeness, since the plan rules it out and the reason should be on
record. The two redefinitions on offer are:

* **"at least two sign changes"** — deletes the only clause that would catch a
  spurious extra crossing. The measured third zero at 1.2204 shows such
  crossings are not hypothetical.
* **"exactly one rising sign change in the design's own gate window [0.10,
  1.0]"** — this is what D-1 tempts you into, and it is arithmetically true of
  both curves (the digitized column's first zero, 0.0656, is itself below 0.10).
  But it *removes the low-x zero from the gate entirely* rather than widening a
  window, and it converts G3a from "the shape matches" to "one feature
  matches". Three documents already say the low-x zero is a poor discriminator
  (`phase_D_gate.md:56-60`, `b1_nuclear.hpp:24-28`, `USAGE.md:271-275`); that is
  an argument for not over-weighting it, not for deleting it.

Widening keeps the clause's content and removes an artefact. That is why it is
defensible and the other two are not. **Whichever is chosen, `design_D_b1_li6.md`
§5.4 must be edited first and the test moved to match** — `PLAN.md:96-97` —
because as things stand the code is a faithful implementation of a
self-contradictory design, and editing only the test would make the code
*unfaithful* to it.

---

## 3. `validation/b1_li6_table.py` — adding a row, what it writes, what it costs

### 3.1 How a new structure-function row is added

Everything a gate row needs is behind one function. `gate_report(unpol=None,
finite_q=True, label="", n_x=N_X)` at **:100-149**:

* **:120-124** builds `lg.DeuteronConvolutionB1.Options()`, sets
  `o.unpol = unpol` **only when `unpol is not None`** (so `None` means the
  library `ToyF2`, matching `Options::unpol`'s "null => ToyF2",
  `b1_nuclear.hpp:555`), sets `o.finite_q_delta = finite_q`, and constructs the
  object. A new structure function therefore needs nothing but an object
  satisfying the `UnpolSF` binding.
* **:126-131** memoizes `x·b₁` per x — the reason the row costs one evaluation
  per grid point and not three.
* **:133** the grid: `sorted(set(x_grid(0.01, 1.0, n_x)) | set(x_grid(0.10, 0.80, n_x)))`.
* **:134-136** `lg.b1_landmarks(xb1, grid)`, `lg.b1_landmarks_of_table()`, and
  the G3b peak over `x_grid(0.10, 0.80, n_x)`.
* **:137-149** returns a dict: `label`, `finite_q`, `zeros`, `zero_slope`,
  `min`, `max`, `integral`, `peak_window`, `g3b_ratio`, and the digitized
  `table` block.

`checklist_item4(n_x=N_X)` at **:152-169** is the list of rows:

```python
def checklist_item4(n_x=N_X):
    rows = [gate_report(None, True,
                        "ToyF2, kappa = sqrt(1+gamma^2)  [Eq. 21, the DEFAULT]",
                        n_x),
            gate_report(None, False, "ToyF2, kappa = 1  [Eq. 17]", n_x)]
    if not lg.HAVE_LHAPDF:
        return rows
    lg.lhapdf_quiet()
    ct18 = lg.LhapdfSF("CT18NLO")
    rows += [gate_report(ct18, True, "CT18NLO, kappa  [Eq. 21]", n_x),
             gate_report(ct18, False, "CT18NLO, kappa = 1", n_x)]
    return rows
```

**To add MSTW2008 LO (task A2), the edit is exactly this shape and nothing
else:** construct the new `UnpolSF` once, append two `gate_report(...)` calls
with `finite_q=True` and `finite_q=False`, and guard it with whatever
availability flag the new class exposes — the `lg.HAVE_LHAPDF` early-return at
**:163-164** is the pattern to copy, and note it currently returns *before* the
CT18 rows, so an MSTW guard must be independent of it or MSTW rows will vanish
on a build without LHAPDF. `emit()` needs no change at all: it iterates
`for r in reports` at **:233-237** and formats whatever it is given, and the
"has no LHAPDF" notice at **:228-230** is the only place that assumes what the
list contains.

The label string is the row's identity in the emitted table (**:234-235**) and
is not parsed anywhere; make it say the δ-function, because
`phase_D_gate.md:14-19` is emphatic that a landmark is meaningless without its
grid and its options.

### 3.2 What the script writes — verified

The module docstring at **:5-7** says:

> Outputs (with `--write`, and it writes NOTHING outside `docs/`):
>     docs/open_items/run_2026-09-02/phase_D_gate.md      -- the gate, section 5
>     docs/open_items/run_2026-09-02/phase_D_numbers.md   -- the numbers, section 7

**That is wrong.** `main()` at **:251-274** is the whole of the write path:

```python
    args = ap.parse_args(argv)
    emit(sys.stdout, args.n_x)
    if args.write:
        # The two documents are hand-written prose around these numbers, so
        # this appends a machine-regenerated appendix rather than overwriting
        # an author's argument with a table.
        out = DOCS / "phase_D_regenerated.md"
        with open(out, "w") as fh:
            fh.write("# Phase D -- regenerated by validation/b1_li6_table.py\n\n" ...)
            emit(fh, args.n_x)
        print("\nwrote %s" % out, file=sys.stderr)
```

`DOCS` is `docs/open_items/run_2026-09-02` (**:41**). So:

* **without `--write`** the script writes **no file at all**; everything goes to
  stdout (**:262**). This is how it was run for the timing below.
* **with `--write`** it writes **exactly one** file,
  `docs/open_items/run_2026-09-02/phase_D_regenerated.md`, truncating it
  (`open(..., "w")`). It never touches `phase_D_gate.md` or
  `phase_D_numbers.md`; the comment at **:264-266** says that is deliberate,
  and the docstring simply was not updated to match.
* `phase_D_regenerated.md` **does not exist in the tree today**, so `--write`
  has not been run here (or its output was not kept). Both prose documents are
  hand-maintained and must be edited by hand after a rerun.

### 3.3 Two more stale numbers in the same file

* **:105-106** `gate_report`'s docstring says the grid is *"(199 points at the
  default --n-x 100…)"*. It is **200**. The two `x_grid` calls share no
  floating-point value, so the `set` union is 100 + 100. The two candidates for
  a collision differ in the last bit: `x_grid(0.01, 1.0, 100)[9]` is
  0.099999999999999991673 against `x_grid(0.10, 0.80, 100)[0]` =
  0.10000000000000000555, and at the other end 0.79999999999999993339 against
  0.80000000000000004441. The script's own emitted header computes the
  length at run time (**:221-226**) and printed **200** in today's run, and
  `phase_D_gate.md:16` says 200.
* **:257-260** `--n-x`'s help says the finite-q rows cost *"~0.14 s per
  point"*. Measured today: **0.242 s per point** (§3.4), which is what the
  `gate_report` docstring at **:110-114**, `b1_nuclear.hpp:573-574` and
  `CONVENTIONS.md:158` all say. The `--help` text is the outlier.

### 3.4 How long a full run takes — measured

Machine: this one, 8 logical CPUs, `source env.sh`, build already in place.

| run | wall | user | notes |
|---|---|---|---|
| `python3 validation/b1_li6_table.py` (defaults, LHAPDF present, 4 gate rows) | **1:39.89** | 99.63 s | 100 % CPU, i.e. single-threaded; 53.5 MB peak RSS; exit 0 |
| the same with `--n-x 2` (4 grid points ⇒ 8 finite-q evaluations) | 4.95 s | 5.14 s | isolates the fixed cost |
| `./build/lipolgen_tests -tc='b1_nuclear*'` | 27.57 s | 93.57 s | 341 % CPU (the layer-3 scan is threaded, `test:104-120`); 19 cases, 3547 assertions, all pass |

Decomposition of the 99.9 s, from those two rows: fixed cost (import, §7 four
terms, the band and knob rows, the densities, the sign gate) ≈ **3.0 s**; the
gate scan is the other ≈ **96.9 s**, which is 2 finite-q rows × 200 grid points
× **0.242 s** per point. The two κ = 1 rows contribute ≈ 0.3 s in total. This
is a direct confirmation of the 0.24 s/point figure and of the factor-300 cost
ratio quoted at `b1_li6_table.py:110-114`.

**Therefore, for planning A2:**

* each additional **finite-q** row at the default `--n-x 100` costs
  **≈ 48 s**; each additional **κ = 1** row costs **≈ 0.15 s**.
* two MSTW rows (one of each) ⇒ ≈ **+48 s**, i.e. a full run of ≈ **2 min 28 s**.
* lowering the *script's* grid floor to 0.001 does not change the point count
  (`x_grid` takes `n`, not a step), so it is free; lowering the *doctest's*
  floor is likewise free, because `linspace(0.001, 1.59, 300)` is still 300
  points.
* today's run reproduced every entry of `phase_D_gate.md:299-305` — G3b ratios
  **0.4400 / 0.2717 / 0.7193 / 0.4532**, second zeros **0.37737 / 0.39230 /
  0.43877 / 0.44920**, minima and maxima to the printed digits. (The script
  prints *all* zeros, so it additionally shows the first zero at 0.02228 /
  0.02231 for the two ToyF2 rows and none at all for the two CT18 rows — the
  fact §2 is about.) The document and the code are in sync as of today.

---

## 4. Every file that carries the b₁ publication ban

`PLAN.md:116-120` says **fourteen** and lists them. The list is short by two:
**`docs/DEVELOPMENT_PLAN.md`** and **`python/bindings.cpp`** also carry the
statement and are not on it.

> **CORRECTED 2026-09-04, after the lift.** This survey was itself short by
> two SITES, and both were missed by every ban-lift pass because of it, so
> they still asserted a failing gate a day after it closed:
> **`python/lipolgen/__init__.py`** (the Sphinx docstring of the public
> `B1_MODELS` constant, which SHIPS in the installed package and contradicted
> `python/bindings.cpp`'s docstring for the same enum) and a **second site in
> `docs/OPEN_ITEMS_SOLUTIONS.md`** (§"8-10. Theory notes", whose pointer to
> §10 landed the reader on the opposite claim). Both are in the table below,
> flagged `MISSED`. The corrected totals are **seventeen files, 48 sites**;
> `OPEN_ITEMS_SOLUTIONS.md` carries it in **8** places, not 7. The lesson the
> table itself teaches: a grep for the *ban* text finds the ban, and a grep
> for "NOT passed" / "fails" finds the FAIL class — this survey ran the first
> and under-ran the second, and the two `__init__.py` / theory-note sites are
> pure FAIL sites with no ban wording in them at all.

**Seventeen files, 48 sites** — the table below is the enumeration and the
counts are its rows. Seven files carry it in three places or more:
`OPEN_ITEMS_SOLUTIONS.md` (8), `phase_D_gate.md` (7),
`USAGE.md` (6), and `design_D_b1_li6.md`, `b1_nuclear.hpp`, `bindings.cpp`,
`test_b1_nuclear.cpp` (3 each). `PHYSICS_CHANNELS.md` carries it exactly twice,
as `PLAN.md:119-120` says.

Legend: **BAN** = an explicit "no ⁶Li number may be published / ships" clause;
**FAIL** = a statement that the gate is not passed, or the 2.27 / 3.68 deficit
quoted as current fact. Both classes have to move when the gate closes; only
the first is the ban proper.

| # | file | anchors | class | what is there |
|---|---|---|---|---|
| 1 | `README.md` | `:49` (as of 66dcda2) | FAIL | "A = 2 magnitude gate is still open, see below" |
| | | `:188-190` | BAN + FAIL | "gate fails its magnitude clause by a factor 2.27–3.68 … no ⁶Li number from the backend may be published until it closes" |
| | | `:62` (as of 66dcda2) | — | the `r1998` follow-up, mentions the gate but is not a ban |
| 2 | `docs/USAGE.md` | `:259` (as of 66dcda2) | FAIL | section heading "⚠ The A = 2 validation gate is NOT fully passed" |
| | | `:262` (as of 66dcda2) | FAIL | "fails its magnitude clause" |
| | | `:295-297` | FAIL | the clause table, G3b row: ratio 0.440, factor 2.27 |
| | | `:271-275` (as of 66dcda2) | FAIL | "G3a's margin, stated" — the [0.02, 1.0] window and the scan floor |
| | | `:284-289` (as of 66dcda2) | FAIL | the attributed residual, CT18 ×1.67 → 0.719 |
| | | `:291-292` (as of 66dcda2) | **BAN** | "no ⁶Li number from this backend may be published while that stands. Open item 10 stays open. The CLI prints the warning on every run." |
| 3 | `docs/PHYSICS_CHANNELS.md` | `:161` (as of adec442) | **BAN** + FAIL | the whole row *"The A = 2 validation gate of that convolution — not passed, and what that forbids"*; also carries seven `b1_nuclear.hpp:NNN` line anchors (see §5) |
| | | `:540` (as of 66dcda2) | **BAN** + FAIL | "fails its own A = 2 magnitude gate by a factor 2.27, so no number from it may be published (§3)" |
| | | | | **this is the document `PLAN.md:119-120` means by "carries it twice"** — confirmed, exactly two |
| 4 | `docs/OPEN_ITEMS_SOLUTIONS.md` | `:215` (as of adec442) | FAIL (**MISSED**) | §"8-10. Theory notes": "The A = 2 validation was done and it **fails on magnitude** — see §10." Anchor as of 2026-09-04; not found by the original survey, and still stale a day after the lift, so a reader following its own pointer to §10 met the opposite claim |
| | | `:19` (as of a98f0a0) | **BAN** + FAIL | item-10 table row: "the A = 2 validation gate FAILS its magnitude clause … so no ⁶Li number from it may be published and the item does not close" |
| | | `:362` (as of 66dcda2) | FAIL | §10 heading "…and the gate is NOT passed" |
| | | `:1051-1065` | **BAN** + FAIL | the "READ THIS BEFORE QUOTING ANY NUMBER BELOW" block |
| | | `:422-431` | FAIL | "The A = 2 gate, measured" clause table (G3a margin at `:413` (as of 66dcda2), G3b FAIL at `:414` (as of 66dcda2)) |
| | | `:1082-1094` | FAIL | the residual budget; `:428` (as of 66dcda2) "the gate's default today … outside G3b" |
| | | `:444` (as of 66dcda2) | **BAN** | heading "The ⁶Li numbers — RECORDED, NOT PUBLISHED (see the warning above)" |
| | | `:556-572` (as of a98f0a0) | **BAN** | "What has to happen before item 10 can close", 7 conditions; `:572` (as of 66dcda2) "Only then re-open G3b. Until it passes, no ⁶Li number ships." |
| 5 | `docs/CONVENTIONS.md` | `:132-133` (as of b1071b1) | FAIL | "quoting the gate at κ = 1 overstated its deficit as a factor 3.68 instead of 2.27" — **no ban**, but the number is stated as current |
| 6 | `docs/DEVELOPMENT_PLAN.md` | `:148` (as of 66dcda2) | FAIL | "ships opt-in behind a magnitude-gate warning" |
| | | `:308-310` | **BAN** | "the b₁(⁶Li) convolution still carries a magnitude-gate warning and publishes no ⁶Li number until the gate closes" |
| | | | | **NOT on `PLAN.md:116-120`'s list of fourteen** |
| 7 | `docs/open_items/run_2026-09-02/design_D_b1_li6.md` | `:1197` | **BAN** | "Nothing about ⁶Li may be quoted, plotted or merged until this passes." |
| | | `:1355-1360` | **BAN** | the *Escalation* clause itself — the source of every other ban line |
| | | `:1624-1629` | **BAN** | Q3: "If it does not close, item 10 stays open and no ⁶Li number ships." |
| 8 | `docs/open_items/run_2026-09-02/phase_D_gate.md` | `:31` | FAIL | the 2026-09-03 revision note, "3.68 rather than 2.27" |
| | | `:38-52` | **BAN** + FAIL | the Verdict table (`:46` G3b FAIL) and `:49-52` "no ⁶Li number may be published while that stands. Item 10 stays open." |
| | | `:54-67` | FAIL | "G3a passes, but read how" and "What closes G3b" |
| | | `:183-198` | FAIL | the two fragile clauses; `:195-198` the G3b measurement |
| | | `:224-244` | FAIL | checklist item 0, "still outside" |
| | | `:299-320` | FAIL | the item-4 table and "One clause regresses" |
| | | `:344-376` | **BAN** + FAIL | residual budget and `:376` "Until it passes, no ⁶Li number ships." |
| 9 | `docs/open_items/run_2026-09-02/phase_D_numbers.md` | `:9-17` | **BAN** + FAIL | the "READ THE GATE FIRST" block, "none of the ⁶Li numbers below may be published" |
| 10 | `docs/open_items/run_2026-09-02/PLAN.md` | `:160-169` | **BAN** + FAIL | the phase-D OUTCOME block, `:166` "no ⁶Li number may be published: item 10 does not close" |
| 11 | `docs/open_items/run_2026-09-02/STATUS.md` | `:13` | FAIL | "A = 2 gate: shape passes, magnitude factor 2.27 low (open)" |
| 12 | `include/lipolgen/b1_nuclear.hpp` | `:11-22` | **BAN** + FAIL | the file's `\file` WARNING block; `:20-21` (as of 66dcda2) "Under design section 5.4 'Escalation' NO 6Li number from `Li6ConvolutionB1` may be published while that stands" |
| | | `:24-28` (as of 66dcda2) | FAIL | G3a's margin: "counting window starts at x = 0.02 … drops below the scan floor" |
| | | `:564-571` (as of 66dcda2) | FAIL | `finite_q_delta`'s doc comment, "moves the deficit from a factor 3.68 to 2.27" |
| 13 | `python/bindings.cpp` | `:1235-1238` (as of 66dcda2) | FAIL | section comment "WARNING: the gate is NOT fully passed (G3b, magnitude)" |
| | | `:1391-1393` (as of 66dcda2) | FAIL | `Li6ConvolutionB1`'s class docstring, "the A = 2 magnitude gate is NOT passed: band every number" |
| | | `:4063-4067` | FAIL | `B1Model::Li6Convolution`'s enum docstring, same statement |
| | | | | **NOT on `PLAN.md:116-120`'s list of fourteen** |
| 14 | `python/lipolgen/cli.py` | `:178-184` | FAIL | `--b1-model`'s `--help` text: "li6-convolution has NOT passed its A = 2 magnitude gate (a factor 2.27 … 3.68 …)" — this is the "`--help`" half of the Escalation clause |
| | | `:383-390` (as of a7b3d18) | **BAN** + FAIL | the run banner; `:389-390` (as of d3ac125) "See docs/…/phase_D_gate.md; no 6Li number from it may be published." |
| 15 | `python/lipolgen/__init__.py` | `:160-162` | FAIL (**MISSED**) | the Sphinx docstring of the public `B1_MODELS` constant: "its A = 2 magnitude gate is NOT passed -- band every number with --b1-band-scale 0/1/2 and never quote one row alone". Anchor as of 2026-09-04. It **ships in the installed package** and describes the same enum as `python/bindings.cpp`, so after the lift the two shipped descriptions of one object said opposite things. **NOT on `PLAN.md:116-120`'s list, and not on this survey's original sixteen either** |
| 16 | `tests/test_b1_nuclear.cpp` | `:487-494` (as of 66dcda2) | FAIL | G3b's comment block, "This FAILS at the default: the ratio is 0.440, a factor 2.27" |
| | | `:504-506` (as of 66dcda2) | FAIL (**executable**) | `CHECK_CLOSE(ratio, 0.439986, 2e-3); CHECK(ratio < 0.5); // the gate's honest state, pinned` |
| | | `:521-533` (as of 66dcda2) | FAIL (**executable**) | the two fragile-clause comments and `CHECK(z[0] > 0.02)`, `CHECK(z[0] - 0.02 < 0.005)`, `CHECK(|z1 − ref| > 0.05)` |
| 17 | `validation/b1_li6_table.py` | `:21-26` | **BAN** + FAIL | module docstring, "READ THE GATE BEFORE QUOTING ANY 6Li NUMBER … no 6Li number from `Li6ConvolutionB1` may be published while that stands" |
| | | `:243-248` | **BAN** + FAIL (**executable**) | the emitted footer: "THE GATE FAILS and no 6Li number may be published (design section 5.4, Escalation)" — this text is *printed into every report*, so it must move or every regenerated report will contradict the gate |

**Not carrying it, checked** *(and see the 2026-09-04 correction above —
`python/lipolgen/__init__.py` was never checked at all, which is how it was
missed)***:** `src/core/b1_nuclear.cpp` (no ban text at all),
`python/tests/test_b1_model.py` (`:14-21` states the rule in the module
docstring but **asserts nothing about it**, so the pytest suite will not break
when the ban lifts), `docs/open_items/run_2026-09-03/PLAN.md` and `STATUS.md`
(they describe the ban as the thing to lift, which is correct either way),
`validation/reference/b1_default_li6.json` (a grep hit on the digits "2.27"
inside data), `data/vmc/deuteron/fdeut.av18` (likewise).

---

## 5. Consequences the implementer should not have to rediscover

1. **The doctest goes red when the gate goes green.** `test:504`, `:505`,
   `:506` pin the failing ratio, and `:531`, `:533` pin the thin margins. All
   five are expected to fail on a successful rerun. Plan the edit with the
   physics change, not after it.
2. **Both grid sites in the doctest**, `test:444` and `test:582`, use
   `linspace(0.01, 1.59, 300)`. Changing one and not the other splits the two
   test cases onto different grids and makes `phase_D_gate.md`'s κ-column
   comparison no longer like-for-like.
3. **`PHYSICS_CHANNELS.md:184-186` carries seven `include/lipolgen/b1_nuclear.hpp:NNN`
   line anchors** (`:551`, `:575`, `:456`, `:603`, `:608`, `:610`, `:629`), and
   `validation/check_physics_channels_links.py` gates them to ±2 lines
   (`:17-18`). Editing the header WARNING block at `b1_nuclear.hpp:11-28`
   shifts every one of them. `run_2026-09-03/PLAN.md:37` already schedules
   `--fix` at the end of each phase; this is the concrete reason it is needed
   here.
   *(Superseded 2026-09-04: the ±2 window is gone. The gate is STRICT by
   default — the anchor must land on the line that declares the name — and
   `b1_nuclear.hpp:11-28` is itself now a fingerprinted RANGE citation, so
   editing that block fails the gate until the row is re-read and
   `--record-ranges` is run. See `phase_C_numbers.md` §§C7–C8.)*
4. **`validation/b1_li6_table.py`'s emitted footer (`:243-248`) is generated
   text.** It is not prose in a document; every future report carries it. It is
   the one ban site that will silently re-appear after the lift if only the
   documents are edited.
5. **`phase_D_regenerated.md` does not exist.** Anyone reconciling the rerun
   against the prose must produce it (`python3 validation/b1_li6_table.py
   --write`) or diff against stdout; it is not in the tree to diff against.
6. **The design document must be edited before the test.** `PLAN.md:96-97`
   says so and it is right for a reason worth stating: today the code is a
   *faithful* implementation of a self-contradictory design (§1.3 D-1, D-2).
   Moving `test:467` without moving `design_D_b1_li6.md:1285-1286` would make
   the code unfaithful to the design and leave the contradiction in place for
   the next reader to rediscover. The same edit should settle D-1 — either
   G3a's landmarks are exempt from the "gate nothing outside [0.10, 0.80]"
   sentence at `:1280-1281`, or that sentence means only G3b. Say which.
7. **Do not let G3c's ∫b₁ dx follow the widened grid** (§2.6a). It is a
   trapezoid of b₁ = (x·b₁)/x; the extension to x = 0.001 moves it +23 % on
   1/x weighting alone, and the design requires it be taken "over the same
   range" as `close_kumano_integral(true)` (`:1300-1301`), which is the
   digitized table's [0.0100, 1.590].
8. **One number needs re-pinning and only one**: `test:515`, `got.x_min`,
   0.237224 → **0.234833** (§2.6b). Everything else in both cases stays inside
   its existing tolerance.
