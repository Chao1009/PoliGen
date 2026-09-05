# Phase A — the A = 2 gate rerun with CDKS's own PDF, every number

Task A2 of `docs/open_items/run_2026-09-03/PLAN.md` §Phase A. Written after the
rerun, from output that is quoted verbatim below with the command that produced
it. Nothing in this file is an estimate; where a number is *not* measured the
sentence says so.

Companion documents: `phase_A_gate_mechanics.md` (the read-only survey that
argued the window change), `../run-2026-09-02/design_D_b1_li6.md` §5.4 (the
criterion, amended today), `phase_D_gate.md` (the previous verdict, now
superseded on G3b and **not yet edited** — see §7), `phase_A_cdbonn.md` (task
A5's offline derivation of the CD-Bonn coefficients).

**§8 was added later, by task A5** (the real CD-Bonn wave function). It does
not change any number in §§0–7 — the gate's default is unchanged — but it
supersedes the `item 5` row of §5's residual budget and closes the residual
that §5 leaves unattributed.

> **STATUS NOTE, added by task A4/A6/A7 on the same day.** Where §0, §7.1 and
> §8.0 item 3 say the ⁶Li publication ban is untouched, that was true when they
> were written and is **no longer true**: the ban was lifted on the strength of
> these measurements, and the **17 files** that carried it now carry the
> verdict and its four conditions instead — 17 is the measured count
> (`phase_A_gate_mechanics.md` §4 enumerates them, and
> `git grep -il 'publication ban'` returns the same 17 files at `ac22331` and
> at `HEAD`; re-measured 2026-09-05).  "Sixteen" here was the mid-task count,
> before the survey finished. **No number in this document
> changed** — only the sentences about what had and had not yet been reconciled
> in the prose. The decision, its conditions and the two author decisions that
> came with it are `docs/OPEN_ITEMS_SOLUTIONS.md` §10 and
> `docs/open_items/run_2026-09-03/STATUS.md`.

---

## 0. The headline

**G3b PASSES.** With MSTW2008 LO — the PDF CDKS themselves used — at CDKS
Eq. (21)'s finite-|q⃗| δ-function, which is `DeuteronConvolutionB1`'s default:

> **max |x·b₁| over [0.10, 0.80] = 9.15096 × 10⁻⁴ against the digitized
> 1.08521 × 10⁻³ → ratio 0.843243, a factor 1.186 low.**
> G3b asks for the ratio in [0.5, 2]. **0.843 is inside, clearing the lower
> edge by a factor 1.69.**

The three PDF families, all at the same δ-function, same grid, same everything
else:

| nucleon PDF | G3b ratio | factor low | G3b verdict |
|---|---|---|---|
| library `ToyF2` | 0.440001 | 2.27 | **outside** [0.5, 2] |
| `LhapdfSF("CT18NLO")` — the stand-in | 0.719432 | 1.39 | inside |
| **`MstwSF()` — MSTW2008 LO, CDKS's own** | **0.843243** | **1.186** | **inside** |

**G3a also passes on MSTW, and passes better than on either stand-in** (§2).
**G3c is reported, not enforced, and is unchanged in character**: still a factor
2.0 below the digitized Close–Kumano integral (§4).

Two things this does **not** say. It does not lift the ⁶Li publication ban:
that is a separate decision with seven conditions listed in
`docs/OPEN_ITEMS_SOLUTIONS.md` §10, and G3b is one of them. And it does not
close item 10; §7 lists what is still un-edited and what is still open.

---

## 1. What was run

Environment: `source env.sh`, build
already in place, 8 logical CPUs. Every command below was run in this tree
today (2026-09-03).

| # | command | wall | what it produced |
|---|---|---|---|
| C1 | `python3 validation/b1_li6_table.py` | 148.8 s | the six-row gate table of §2.2 |
| C2 | `./build/lipolgen_tests -tc='b1_nuclear*'` | 49.5 s | 20 cases, 5686 assertions, all pass |
| C3 | `./build/lipolgen_tests` | 99.3 s | **368 cases / 17206076 assertions / 1 skipped, 0 failed** |
| C4 | `python -m pytest python/tests -q` | 34.3 s | **236 passed** |
| C5 | a multiprocessing replica of the doctest's layer-3 scan (scratch script), 10 rows × 300 points | 64.2 s | §2.1's both-grid table, and the values now pinned in the doctest |
| C6 | a short `f1_nucleon` / `f2p` probe (scratch) | 4 s | §5's input-ratio table |

C1's full output, including the four terms, the band rows, the densities and
the α–d sign gate, is reproducible by rerunning it; only the layer-3 section is
transcribed here.

---

## 2. G3a — shape

### 2.1 The doctest's grid, `linspace(0.001, 1.59, 300)` — the gated numbers

These are the numbers `tests/test_b1_nuclear.cpp` now pins. Produced by C5 and
confirmed assertion-by-assertion by C2.

| | ToyF2, κ | ToyF2, κ = 1 | CT18NLO, κ | **MSTW, κ** | MSTW, κ = 1 |
|---|---|---|---|---|---|
| all sign changes found | 0.022040(−), 0.377374(+), 1.220437(−) | 0.022060(−), 0.392320(+), 1.142281(−) | 0.009795(−), 0.438761(+), 1.197722(−) | **0.027886(−), 0.495239(+), 1.217660(−)** | 0.027946(−), 0.507853(+), 1.137076(−) |
| count in **(0, 1.0]** — G3a's clause | **2 ✓** | 2 ✓ | **2 ✓** | **2 ✓** | 2 ✓ |
| count in [0.02, 1.0] — the *old* clause | 2 ✓ | 2 ✓ | **1 → ABORT** | **2 ✓** | 2 ✓ |
| first zero **falling** | ✓ | ✓ | ✓ | ✓ | ✓ |
| \|z₀ − 0.065645\|, tol **0.08** | 0.043605 ✓ | 0.043585 ✓ | 0.055850 ✓ | **0.037759 ✓** | 0.037699 ✓ |
| second zero **rising** | ✓ | ✓ | ✓ | ✓ | ✓ |
| \|z₁ − 0.457177\|, tol **0.10** | 0.079803 ✓ | 0.064857 ✓ | 0.018416 ✓ | **0.038062 ✓** | 0.050676 ✓ |
| peak position in [0.5, 1.0] | 0.755642 | 0.739699 | 0.734385 | **0.771585** | 0.755642 |
| \|x_peak − 0.765663\|, tol **0.10** | 0.010021 ✓ | 0.025964 ✓ | 0.031278 ✓ | **0.005922 ✓** | 0.010021 ✓ |
| **G3a** | **PASS** | PASS | **PASS** (fails the old clause) | **PASS** | PASS |

**CDKS's own PDF reproduces CDKS's own figure better than either stand-in on
every G3a landmark.** The second zero misses by 0.038 against ToyF2's 0.080,
and the peak position by 0.0059 against ToyF2's 0.0100 and CT18NLO's 0.0313.
That is a non-trivial consistency check: nothing was tuned, the only thing that
changed is which published grid the nucleon F₁ was read from.

The reference column (`b1_landmarks_of_table()`, computed not typed):
zeros **0.065645** (falling), **0.457177** (rising); minimum **−1.768140e−4**
at **0.332355**; maximum **+1.085210e−3** at **0.765663**; ∫b₁ dx
**+4.592003e−4**.

### 2.2 The script's grid — C1's table verbatim

Union of `x_grid(0.001, 1.0, 100)` and `x_grid(0.10, 0.80, 100)`, 200 points.
It is a **different grid** from the doctest's (coarser, and it stops at 1.0 not
1.59), so landmarks differ in the 4th digit and the integrals are over a
different range. That is a grid difference and nothing else.

```
| configuration | zeros (G3a) | minimum | maximum | int b1 dx (x in [0.01, 1.0]) | G3b ratio |
|---|---|---|---|---|---|
| ToyF2, kappa = sqrt(1+gamma^2)  [Eq. 21, the DEFAULT] | 0.02214, 0.37735 | -3.3089e-05 @ 0.2343 | +4.7751e-04 @ 0.7576 | +2.0666e-04 | **0.4400** |
| ToyF2, kappa = 1  [Eq. 17] | 0.02216, 0.39233 | -3.1211e-05 @ 0.2343 | +2.9483e-04 @ 0.7376 | +1.0555e-04 | **0.2717** |
| CT18NLO, kappa  [Eq. 21] | 0.01010, 0.43877 | -1.8676e-04 @ 0.2980 | +7.8062e-04 @ 0.7364 | +2.0140e-04 | **0.7193** |
| CT18NLO, kappa = 1 | 0.01010, 0.44918 | -1.7025e-04 @ 0.2980 | +4.9178e-04 @ 0.7152 | +5.3563e-05 | **0.4532** |
| MSTW2008 LO, kappa  [Eq. 21]  <- CDKS's own PDF | 0.03029, 0.49523 | -2.3458e-04 @ 0.3545 | +9.1509e-04 @ 0.7717 | +2.0979e-04 | **0.8432** |
| MSTW2008 LO, kappa = 1 | 0.03035, 0.50786 | -2.0832e-04 @ 0.3545 | +5.6457e-04 @ 0.7505 | +5.4988e-05 | **0.5202** |
| *digitized CDKS (b1_landmarks_of_table)* | 0.06564, 0.45718 | -1.7681e-04 @ 0.3324 | +1.0852e-03 @ 0.7657 | +4.5920e-04 | *1* |
```

The G3b ratios agree with §2.1's to 4 digits on every row (0.4400 / 0.2717 /
0.7193 / 0.4532 / 0.8432 / 0.5202), which is the cross-check that the two
independent scans are measuring the same thing.

### 2.3 How well the low-x zero is actually resolved — a caveat, measured

The same CT18NLO first zero comes out **0.01010** on the script's 200-point
grid, **0.009795** on the doctest's 300-point grid, and **0.009396** on the
dedicated log probe of `phase_A_gate_mechanics.md` §2.3. `b1_landmarks` sees a
zero only *between two consecutive grid points* and places it by linear
interpolation, so on a uniform grid this landmark is resolved to a few per
cent at best. Its ±0.08 tolerance swallows all of that (the widest miss is
0.056), so it changes no verdict — but **the low-x zero should not be quoted
to five digits anywhere.** The same caveat applies to MSTW's 0.0279 / 0.0303.

---

## 3. G3b — magnitude, the clause the phase exists for

`max |x·b₁|` over x ∈ [0.10, 0.80] against the **raw** digitized column's
maximum 1.085210 × 10⁻³. Doctest grid (C2/C5); the script's grid (C1) agrees
to 4 digits.

| row | max\|x·b₁\| | ratio | factor low | in [0.5, 2]? |
|---|---|---|---|---|
| ToyF2, κ = 1 | 2.94839e−4 | 0.271689 | 3.68 | no |
| ToyF2, **κ** (the default) | 4.77494e−4 | 0.440001 | 2.27 | no |
| CT18NLO, κ = 1 | — (script grid: 0.4532) | 0.4532 | 2.21 | no |
| CT18NLO, **κ** | 7.80734e−4 | 0.719432 | 1.39 | yes |
| MSTW, κ = 1 | 5.64669e−4 | **0.520332** | 1.92 | **yes, by 4 %** |
| **MSTW, κ (the gate's default)** | **9.15096e−4** | **0.843243** | **1.186** | **YES** |

**Why the MSTW/κ row is the verdict and not one row among six.** Design §5.4's
Escalation clause is *"If **after the checklist** the peak ratio is still
outside a factor of 2, stop and hand back"*. Checklist item 4 **is** the
nucleon PDF, it is the dominant term of the residual by the design's own
budget, and CDKS built the Fig. 4 curve on MSTW2008 LO. Reading G3b off ToyF2
would be reading it before the checklist; reading it off CT18NLO would be
reading it off a stand-in that this phase exists to remove.

**The honest reading of the κ = 1 row.** At Eq. (17)'s δ-function MSTW gives
0.5203 — inside the window by 4 % of its own value. So a meaningful part of
G3b's pass is carried by the choice of Eq. (21) over Eq. (17). That choice is
defensible and is CDKS's own (Eq. 21 is their exact definition, Eq. 18's "≈"
is what makes Eq. 17 the approximation), and it was made and defended on
2026-09-03 independently of this measurement. But it should be stated in the
same breath as the pass, not buried: **G3b passes comfortably at Eq. (21) and
marginally at Eq. (17).**

---

## 4. G3c — Close–Kumano, reported not enforced

∫b₁ dx from this kernel, on the **x ≥ 0.01 sub-grid** of the doctest's scan
(i.e. [0.011629, 1.59]), against `close_kumano_integral(true)` = **+4.592003e−4**
(the digitized CDKS column over its own [0.0100, 1.590]) and
`close_kumano_integral(false)` = **+5.91474e−3** (Miller).

| row | ∫b₁ dx, x ≥ 0.01 | ratio to digitized | ∫ over the **whole** scan (NOT G3c) |
|---|---|---|---|
| ToyF2, κ | +2.14607e−4 | 0.467 | +2.65316e−4 |
| ToyF2, κ = 1 | +1.07152e−4 | 0.233 | +1.57861e−4 |
| CT18NLO, κ | +2.09620e−4 | 0.456 | +2.47427e−4 |
| **MSTW, κ** | **+2.24896e−4** | **0.490** | +2.66866e−4 |
| MSTW, κ = 1 | +5.79583e−5 | 0.126 | +9.99290e−5 |

Nothing here is enforced (`sf.hpp`: "Reported, not enforced"), and none of the
three is zero. **MSTW moves G3c barely at all** — 0.467 → 0.490 — which is
worth stating plainly: the PDF that fixes the *peak* does not fix the
*integral*. The integral is dominated by the low-x region where all three
curves sit well below the digitized column, and by the sign structure, and the
PDF is not the lever there.

**The whole-scan column is not G3c** and must never be quoted as it. It is
18 % (MSTW, κ) to 72 % (MSTW, κ = 1) larger purely because `integral_b1` is a trapezoid of b₁ = (x·b₁)/x and
the extension to the new scan floor 0.001 adds a region where 1/x ≈ 10³. Design
§5.4's G3c now says this explicitly; both numbers are pinned in the doctest so
they can never be confused again.

---

## 5. The residual budget, restated against 0.843

The design's budget (`design_D_b1_li6.md` §5.4 G3b, and `phase_D_gate.md`
"Residual budget") starts from the κ = 1 ToyF2 column and multiplies. Redone
on the doctest grid with the measured MSTW rows:

| piece | measured effect on the peak | direction |
|---|---|---|
| starting point: ToyF2, κ = 1 | ratio **0.271689** | — |
| item 0 — finite-\|q⃗\| δ-function, κ = 1 → √(1+γ²) | **×1.6195** (ToyF2) / **×1.6206** (MSTW) | closes — **IN** (the default) |
| item 4 — nucleon PDF, ToyF2 → **MSTW2008 LO** | **×1.9152** (at κ = 1) / **×1.9165** (at κ) | closes — **NOW IN** |
| *(for comparison: ToyF2 → CT18NLO)* | ×1.6680 (at κ = 1) | the stand-in, superseded |
| *(CT18NLO → MSTW at κ)* | ×1.1721 | what this phase bought |
| item 1 — target mass (CDKS Eq. 22, already in) | ×1.487 at x = 0.8, κ (1.520 at κ = 1) | in |
| item 2 — R (`r1998`, CDKS's own; already in) | ×0.98 | in |
| item 5 — wave function, AV18 → rescaled to P_D = 4.85 % | ×0.88 | **opens** |
| figure digitization | not quantified | — |
| **the gate's default, MSTW + κ** | **ratio 0.843243 (a factor 1.186 low)** | **INSIDE G3b** |

*Superseded on the `item 5` row by §8*: that ×0.88 is the **rescaling proxy**'s.
The real CD-Bonn wave function is **×1.18631** on this row and **closes** — it
supplies the whole of the 1.186 that the next paragraph calls unattributed.

**The two large items are independent to 0.07 %**, which is the check that
there is no third unexplained factor at the peak:

```
0.271689 x 1.91519 x 1.61951 = 0.842682   predicted
                               0.843243   measured   (0.07 % overlap)
```

**What is left unattributed: a factor 1.186, i.e. the computed peak is 18.6 %
below the digitized one.** The design budgeted four sources for exactly this
residual and two of them remain live and unquantified — AV18 vs CD-Bonn (the
design says ~20 % on P_D and more on the tail; `phase_A_cdbonn.md` is the
separate task), and the fact that the target curve is a *figure
digitization*. An 18.6 % residual is inside that band. **It is not evidence
that the remaining physics is understood; it is the statement that the gate can
no longer distinguish it, which is what the design says a factor 2 is for.**

### 5.1 Why the input ratio is 1.40 but the output ratio is 1.17

The isoscalar F₁ the convolution is actually fed (`f1_nucleon(x, 2.5)`, C6),
and the F2p ratio for reference:

| x | F₁ ToyF2 | F₁ CT18NLO | F₁ MSTW | MSTW/CT18 on F₁ | MSTW/ToyF2 on F₁ | MSTW/CT18 on F2p |
|---|---|---|---|---|---|---|
| 0.10 | 1.388924 | 1.549933 | 1.381344 | 0.89123 | 0.99454 | 0.88493 |
| 0.30 | 2.824588e−1 | 4.209263e−1 | 4.163149e−1 | 0.98904 | 1.47390 | 0.98386 |
| 0.50 | 8.353448e−2 | 1.311235e−1 | 1.517377e−1 | 1.15721 | 1.81647 | 1.16833 |
| 0.70 | 1.694132e−2 | 2.408928e−2 | 3.281082e−2 | 1.36205 | 1.93673 | 1.41759 |
| 0.736 | 1.148925e−2 | 1.585638e−2 | 2.222510e−2 | 1.40165 | 1.93442 | 1.46783 |
| 0.770 | 7.617219e−3 | 1.018464e−2 | 1.465636e−2 | 1.43907 | 1.92411 | 1.51520 |
| 0.80 | 5.053784e−3 | 6.549908e−3 | 9.639765e−3 | 1.47174 | 1.90744 | 1.55631 |

MSTW/ToyF2 on F₁ is **1.91–1.94 for x ≥ 0.7** and the peak x·b₁ ratio moved by
**×1.9165** — the two agree to 1 %, so at ToyF2 the peak is essentially a
pointwise read of F₁. MSTW/CT18NLO on F₁ is **1.40–1.44** at the peak but the
peak x·b₁ ratio moved only **×1.1721**. The convolution is not a pointwise
multiplication — x·b₁(x) integrates F₁(x/y) against the light-cone densities
over y, and the peak itself *moves* (0.7344 → 0.7716) — so the two need not
agree, and the honest statement is that **the gain over CT18NLO is smaller
than the input ratio at the peak x suggests.**

### 5.2 The curve point by point, against the raw digitized column

x·b₁ at the default δ-function (C6), against `2 × CdksB1` (the raw column, the
same normalisation G3b uses):

| x | digitized (raw) | ToyF2 | CT18NLO | **MSTW** | MSTW / digitized |
|---|---|---|---|---|---|
| 0.10 | −1.7289e−5 | −2.0339e−5 | −4.0018e−5 | **−2.2807e−5** | 1.32 (too negative) |
| 0.20 | −8.7797e−5 | −3.2268e−5 | −1.2473e−4 | **−9.8083e−5** | 1.12 |
| 0.30 | −1.6897e−4 | −2.8569e−5 | −1.8679e−4 | **−2.0772e−4** | 1.23 |
| 0.50 | +1.4833e−4 | +1.3260e−4 | +1.9068e−4 | **+1.6465e−5** | **0.11** |
| 0.70 | +9.9670e−4 | +4.5125e−4 | +7.6459e−4 | **+8.2208e−4** | 0.82 |
| 0.80 | +1.0665e−3 | +4.6372e−4 | +7.2253e−4 | **+9.0101e−4** | 0.84 |

Two readings that are not in the headline.

1. **The low-x tail, which `phase_D_gate.md` recorded as the one thing still
   unexplained, largely closes with MSTW.** That document quoted CT18 + κ
   giving x·b₁(0.10) = −3.8e−5 against the digitized −1.73e−5, a factor 2.2.
   With MSTW it is −2.28e−5, a factor **1.32**. The mid-x dip likewise: ToyF2
   was 5.3× too shallow at the minimum (−3.31e−5 against −1.768e−4), MSTW is
   1.33× too **deep** (−2.346e−4). The sign of the miss has flipped, which is
   what "the gate can no longer distinguish these" looks like.
2. **x = 0.50 is the one point that got worse, and badly.** MSTW gives
   +1.65e−5 against +1.48e−4, a factor 9 low, because MSTW's second zero sits
   at 0.4952 and x = 0.50 is 0.005 past it while the digitized curve crossed at
   0.4572 and is already well positive. This is not a magnitude failure, it is
   the zero-position miss of 0.038 read at the worst possible x. G3a gates the
   zero *position* (±0.10, cleared by 0.062) and G3b gates the *peak*, so no
   clause sees it — but a reader comparing curves at x = 0.5 will, and should
   be told why.

---

## 6. What changed in the tree, and what it cost

### 6.1 The design (`design_D_b1_li6.md` §5.4) — amended first, as `PLAN.md:96-97` requires

* **G3a's counting window `[0.02, 1.0] → (0, 1.0]`**, and the **scan floor
  0.01 → 0.001 stated in the clause** instead of being left to whatever grid
  the implementation picked. Ceiling unchanged at 1.0 and load-bearing: every
  configuration measured has a real third crossing at 1.14–1.22.
* **G3c** now says explicitly that "the same range" means the **x ≥ 0.01
  sub-grid**, because `integral_b1` is a trapezoid of (x·b₁)/x.
* **The [0.10, 0.80] "gate window" sentence** is now scoped: it governs **G3b
  only**. G3a is a shape clause whose landmarks (0.0656, 0.766) are themselves
  outside that window, so it could not be stated there; G3c is over the
  digitized range. The design previously said both things and the code
  implemented one.

**Honesty note, and it matters.** The window change was argued in
`phase_A_gate_mechanics.md` from CT18NLO, whose low-x zero at 0.0098 falls
below *both* the old floor and the old window and therefore aborted
`REQUIRE(z.size() == 2)`. **It turned out not to be needed for the measurement
that closes the gate.** MSTW's low-x zero is at **0.0292** on the old grid and
**0.0279** on the new one — above the old floor and above the old counting
window both — so G3a passes on MSTW with or without the amendment (§2.1, row
"count in [0.02, 1.0]"). The amendment was applied anyway, for two reasons
that do not depend on which PDF is used: the inconsistency it removes is real
(a counting window narrower than the position tolerance it brackets makes that
tolerance dead code), and the CT18NLO row is now itself pinned in the doctest,
where it *does* need it. **The step-3 premise of the task — "G3a will abort
before G3b is reached" — was measured on CT18NLO and does not transfer to
MSTW.**

### 6.2 The doctest (`tests/test_b1_nuclear.cpp`)

Grid `linspace(0.01, 1.59, 300) → linspace(0.001, 1.59, 300)` at **both** sites
(the layer-3 case and the checklist 0/1/2 case; changing one and not the other
would split them onto different grids). Counting filter
`>= 0.02 && <= 1.0 → > 0.0 && <= 1.0`. G3c taken on the x ≥ 0.01 sub-grid,
reusing the same evaluations.

Pins that **had** to move (outside their existing tolerance):

| pin | was | now | why |
|---|---|---|---|
| `got.x_min` | 0.2372 | **0.23483** | the argmin hops one grid point on a flat stretch; its *value* moves 3.2e−5 relative |
| `got.integral_b1` (G3c) | 2.15370e−4 | **2.14607e−4** | the sub-grid's first point is 0.011629, not 0.01 — a 0.35 % endpoint relocation, just outside the 3e−3 tolerance |
| `l1.integral_b1` (G3c, κ = 1) | 1.07915e−4 | **1.07152e−4** | same cause, 0.71 % |

Pins refreshed to the new measured value although the old one still held
(`win`, `ratio`, `xb1_min`, `x_max`, `xb1_max`, `z[0]`, `z[1]`, `lr.zeros[1]`,
`l1.zeros[1]`, `lr.xb1_max`, `l1.xb1_max`, `l1.zeros[0]`, `l1.x_max`, the
κ = 1 G3b ratio) — all moved in the 5th–6th digit only, and the file's
convention is a tight pin, so leaving stale values would have hidden the shift.

Assertions **removed**, and why: `CHECK(z[0] > 0.02)` and
`CHECK(z[0] - 0.02 < 0.005)` pinned the margin between the first zero and a
counting-window floor that no longer exists. Leaving them would assert a bound
the design no longer states.

Assertions **added**: the whole-scan `got.integral_b1` and `l1.integral_b1`
(pinned *next to* the G3c sub-grid values so the two can never be confused),
and a whole new case.

**New case: `b1_nuclear T1 gate checklist item 4: the nucleon PDF (G3a/G3b/G3c)`.**
Runs the full layer-3 clause set on MSTW2008 LO at both δ-functions and on
CT18NLO at the default, on the same grid, and pins every landmark plus the
three G3b ratios. It is `#if defined(LIPOLGEN_HAVE_PYTHIA8) || defined(LIPOLGEN_HAVE_LHAPDF)`
and each half additionally guarded, with a runtime check for the MSTW grid
file — the one place in T1 where a skip is right, because these rows need an
*optional tier* while layers 0–3 need only `data/vmc` and must go red without
it. *(**Amended 2026-09-04.** The runtime check was a silent one: with the grid
absent the case still reported "1 passed" and the verdict clause was simply not
evaluated, while every doc in the tree said PASSED in that build too. The MSTW
half now lives in its own case, `b1_nuclear T1v gate checklist item 4, THE
VERDICT ROW: MSTW2008 LO at CDKS Eq. (21) -- G3a/G3b/G3c`, decorated with
`doctest::skip(!mstw_grid_present())` so doctest's own tally counts it SKIPPED
instead of passed — doctest reports a skipped case in the COUNT only, never by
name, so the name and the reason come from `b1_nuclear T1r`, which runs in
**every** build and prints whether the verdict row was measured. This case keeps the CT18NLO half and is now
`#ifdef LIPOLGEN_HAVE_LHAPDF`.)* The scan is threaded with **one structure function per thread**: neither
`MstwSF` (`Pythia8::PDF::xf` is non-const and memoizes the last (x, Q²)) nor
`LhapdfSF` is safe to share, and the Options are built in the calling thread so
the grid reads are serialized. Cost: **+21.6 s** on the suite (77.7 → 99.3 s).

### 6.3 The validation script (`validation/b1_li6_table.py`)

Two MSTW rows added to `checklist_item4()` beside the two ToyF2 and two CT18NLO
rows, guarded on `HAVE_PYTHIA8` + a construction attempt and **independently of
`HAVE_LHAPDF`** (the LHAPDF guard used to be an early `return`, which would
have swallowed the MSTW rows on a build without LHAPDF). Grid floor 0.001,
`integral` on the x ≥ 0.01 sub-grid with the whole-scan value kept beside it.
The emitted table gained an ∫b₁ dx column. The emitted footer — which is
*generated text that goes into every future report* — no longer says "THE GATE
FAILS"; it names the verdict row, prints whether the ratio is inside [0.5, 2],
and states that the ⁶Li publication ban is a separate decision that nothing in
the report lifts. Two stale numbers in the docstrings were corrected on the way
past (the grid is 200 points, not 199).

Cost: **148.8 s** (C1). `phase_A_gate_mechanics.md` §3.4 measured 99.9 s for
the four-row version; that figure was not re-measured here, and the difference
is one extra finite-|q⃗| row at ≈ 48 s, which matches.

### 6.4 Suites

| suite | before this task (with A1 in the tree) | after |
|---|---|---|
| `./build/lipolgen_tests` | 367 cases / 17203937 assertions / 1 skipped | **368 cases / 17206076 assertions / 1 skipped, 0 failed** |
| `python -m pytest python/tests -q` | 236 passed | **236 passed** |

(The brief's stated baseline of 363 / 17203863 / 1 and 211 passed predates task
A1, which added `tests/test_mstw_sf.cpp` — 4 cases — and
`python/tests/test_mstw_sf.py` — 25 tests. 363 + 4 = 367 and 211 + 25 = 236,
so the two baselines reconcile exactly.)

`b1_nuclear T9`, the 600-row rtol-1e−12 gate against
`validation/reference/*.json`, passes unchanged: **no library behaviour was
touched by this task**, only test grids and reporting.

---

## 7. What is NOT done, and what a reader must not conclude

1. **The ⁶Li publication ban is untouched.** **17 files** carry it or a
   gate-failure statement — `phase_A_gate_mechanics.md` §4 enumerates the 46
   sites in them, and `git grep -il 'publication ban'` returns the same 17 at
   `ac22331` and at `HEAD` (re-measured 2026-09-05; this line said "Sixteen",
   the count as of the middle of the survey). This task edited exactly two of them, and only where leaving them
   would have made the *code and its own generated output* contradict the
   measurement: `validation/b1_li6_table.py`'s module docstring and its
   emitted footer. **`phase_D_gate.md`, `OPEN_ITEMS_SOLUTIONS.md` §10,
   `USAGE.md`, `README.md`, `PHYSICS_CHANNELS.md`, `b1_nuclear.hpp`,
   `bindings.cpp`, `cli.py` and the rest still say the gate fails.** They are
   now wrong on G3b and right on everything else. Reconciling them is a
   separate, larger edit with a policy decision inside it.
2. **G3b passing is one of the seven conditions** listed at
   `OPEN_ITEMS_SOLUTIONS.md` §10 ("The seven close conditions, and where each
   one landed") for item 10 to close. The others —
   CD-Bonn, the `r1998` default question, the VMC 5 % spread — are not
   addressed here. (**CD-Bonn was addressed afterwards, by task A5: see §8.**
   The `r1998` question and the VMC spread remain open.)
3. **The κ = 1 column clears G3b by 4 %.** Any statement of the form "the gate
   passes" that does not also say "at Eq. (21)'s δ-function" is overselling it.
4. **G3c did not improve** (0.467 → 0.490 of the digitized integral) and is not
   enforced. The peak agreeing does not make the integral agree.
5. **x = 0.50 got worse**, by a factor 9 low, for the zero-position reason of
   §5.2. No clause of G3 sees it.
6. **The low-x zero is resolved to a few per cent at best** on any of these
   grids (§2.3). Do not quote 0.0279 or 0.0098 to more than two significant
   figures.
7. **CT18NLO at κ = 1 on the doctest grid was not measured** — only on the
   script's grid (0.4532). It is the one cell of §3's table taken from a
   different grid, and it is labelled as such.
8. **`phase_D_regenerated.md` still does not exist**; `--write` was not run.
   The C1 output above is the transcript to diff against.

---

## 8. Task A5 — the real CD-Bonn wave function in the gate

Added after §§0–7 were written, by task A5 of the same `PLAN.md`. The
coefficients and the offline verification are `phase_A_cdbonn.md`; **this
section is only what the repository's own gate measures now that they are
wired in.** Everything below was produced by `build/lipolgen_tests`
(checklist item 5, `tests/test_b1_nuclear.cpp`) or by a probe using the same
`gate_row` code on the same grid — `linspace(0.001, 1.59, 300)`, Q² = 2.5,
CDKS Eq. (21)'s δ-function unless a row says κ = 1.

**The default did not move.** `DeuteronConvolutionB1::Options::wave` defaults
to `kFdeutFile`, so §§0–7's rows, the item-4 pins and every published gate
number keep their meaning exactly. CD-Bonn is `wave = kCdBonn`.

### 8.0 The headline

**With CD-Bonn *and* MSTW2008 LO — CDKS's own wave function and CDKS's own
PDF — the kernel reproduces the digitized CDKS Fig. 4 on every G3 landmark.**

| landmark | digitized CDKS | AV18 (§0's verdict row) | **CD-Bonn** |
|---|---|---|---|
| **G3b peak ratio** | 1 | 0.843243 | **1.000338** |
| peak position x | 0.765663 | 0.771585 | **0.766271** |
| first zero | 0.065645 | 0.027886 | **0.064129** |
| second zero | 0.457177 | 0.495239 | **0.457018** |
| dip x·b₁ | −1.76814e−4 | −2.34606e−4 | **−1.769065e−4** |
| dip position | 0.33236 | 0.35706 | **0.33049** |
| **G3c** ∫b₁dx (x ≥ 0.01) | 4.59200e−4 | 2.24896e−4 | **4.48580e−4** |

Nothing was tuned. The coefficients are Machleidt's Table XX typed once in
`cluster.hpp`; the PDF is the MSTW grid PYTHIA ships; the kernel is byte-identical
to the one §0 was measured with.

**Three things this does not say, and they matter more than the table.**

1. **A 0.03 % agreement is better than the reference deserves.** The target is
   a *digitization of a published figure*. §2.3 already measured that this
   gate cannot resolve the low-x zero better than a few per cent. Read the
   row as "the residual is now below the digitization error", never as
   three-digit agreement with CDKS.
2. **It is specific to CD-Bonn + MSTW.** On CT18NLO the same swap improves the
   peak and the low-x zero and *degrades* the dip and the second zero (§8.3).
   The clean picture appears only when both of CDKS's inputs are used at once,
   which is what one would expect and is also why it is not a coincidence —
   but a reader must not generalise it to "CD-Bonn is better".
3. **It does not lift the ⁶Li publication ban** and does not close item 10.
   G3b was already passing (§0); this makes the residual smaller, not the
   remaining conditions fewer.

### 8.1 The wave function itself, measured through the repository

`cdbonn_wave()` against Machleidt's own Table XV (`tab_deu`). The momentum
moments are **closed form** — with (2/π)∫dp p²/[(p²+a²)(p²+b²)] = 1/(a+b) the
norms are the double sums Σᵢⱼ CᵢCⱼ/(mᵢ+mⱼ) — so there is no quadrature error
anywhere in this table.

| quantity | computed here | CD-Bonn Table XV | gap |
|---|---|---|---|
| (2/π)∫dp p²(u²+w²) | 0.99999982615384 | 1 | 1.74e−7 |
| **P_D** | 0.048562073644234 → **4.8562 %** | 4.85 % | 6.2e−5 |
| **η = A_D/A_S = D₁/C₁** | 0.0255713786530431 | 0.0256(4) | 2.9e−5 |
| **A_S = C₁** | 0.88472985 | 0.8846(9) | 1.3e−4 |
| **Q_d** (`alpha_d_quadrupole_fm2`) | **+0.270178 fm²** | 0.270 fm² | 1.8e−4 |

The gaps in P_D, η and A_S are **the fit residual, not an error**: this
parameterisation is a fit to Machleidt's numerical wave function (he quotes
{∫dr[u−u_a]²}^½ = 2.2e−4, {∫dr[w−w_a]²}^½ = 1.1e−4), so 4.8562 % genuinely does
not round to his 4.85 %. Saying that is why the doctest carries two clauses
per quantity: a loose one against the *published* value and a tight pin (1e−9
relative) against the *computed* one, the pin being what actually catches a
mistyped digit.

**How sharp the pin is, measured.** Perturbing each of the eighteen published
coefficients by one unit in its last printed figure moves the norm by, in units
of 1e−9 relative: C₁ 25, C₂ 4.6, **C₃ 0.21**, C₄ 79, C₅ 42, C₆ 243, C₇ 143,
C₈ 83, C₉ 45, C₁₀ 18, D₁ 16, D₂ 4.8, **D₃ 0.96**, D₄ 23, D₅ 70, D₆ 21, D₇ 55,
D₈ 9.8. At the 1e−9 pin a single wrong digit anywhere is caught with **one
honest exception**: the last figure of C₃ = −0.44114404e−01, an absolute 1e−9
on the smallest coefficient in the table, is worth 2.1e−10 and is not
resolvable in double precision by any observable. A slip in its *seventh*
figure is caught (2.1e−9).

**The derived coefficients**, computed from the r → 0 boundary conditions and
never typed (Eqs. (D23)/(D24) in the equivalent sum-rule form): C₁₁ =
42.260718143999839, D₉ = 1374.064272638856, D₁₀ = −630.43718361517006,
D₁₁ = 120.68764280431213, reproducing `phase_A_cdbonn.md` §3.1's independently
obtained values to 1.4e−14 relative. The four constraint residuals are
0, −3.46e−14, −1.88e−12, −1.18e−10 against sums of scale 1e0/1e3/1e5 — i.e.
1.7e−15 relative, inside Machleidt's "about 15 decimal digits" demand.

**The sign**, which is the `(−i)^L` trap `b1_nuclear.hpp` warns about and which
CD-Bonn walks into (his Eq. (D13) is printed with a bare j_L kernel for both L,
which is inconsistent at L = 2 with his own (D20)/(D22)). Settled by the
repository's own gate, all three from one code path:

| (u, w) | `alpha_d_quadrupole_fm2` | reference |
|---|---|---|
| `fdeut.av18` (control) | **+0.269362 fm²** | its own header `qm` = 0.269673 |
| CD-Bonn, w = −ψ₂ᵃ (what is implemented) | **+0.270178 fm²** | Table XV, Q_d = 0.270 |
| CD-Bonn, w = +ψ₂ᵃ | **−0.303719 fm²** | wrong sign *and* wrong magnitude |

**Against AV18, on the shared 0.1 fm⁻¹ grid** — the two are directly
comparable, both normalised to 1 in the same convention:

| | AV18 | CD-Bonn |
|---|---|---|
| P_D (trapezoid on the grid) | 0.0575985 (header 0.057599) | 0.0485621 |
| ∫k²(u²+w²)dk | 0.9999764 | 0.9999780 |
| w(CD-Bonn)/w(AV18) at 0.1 / 1.0 / 5.0 fm⁻¹ | — | **1.0203 / 0.9888 / 0.3597** |
| first S node | between 2.0 and 2.1 fm⁻¹ | **between 2.3 and 2.4** |

The third row is the point, and it is what condemns the item-5 **rescaling
proxy**: the proxy is *one number* multiplying w everywhere, and the real ratio
runs from 1.02 to 0.36. CD-Bonn is not a softer AV18; it is larger where the D
wave is largest and much smaller in the tail.

### 8.2 The G3 clauses, AV18 vs CD-Bonn, all rows

Every row on the doctest grid at Q² = 2.5. "nz" is the number of sign changes
in G3a's (0, 1.0] counting window.

| PDF | κ | wave | nz | zeros | dip | peak | **G3b** | G3c |
|---|---|---|---|---|---|---|---|---|
| ToyF2 | Eq. 21 | AV18 | 2 | 0.022040, 0.377374 | −3.30902e−5 @ 0.23483 | 4.774939e−4 @ 0.755642 | 0.440001 | 2.14607e−4 |
| ToyF2 | Eq. 21 | **CD-Bonn** | **0** | *none* | — | 5.737089e−4 @ 0.750328 | **0.528662** | 4.49904e−4 |
| CT18NLO | Eq. 21 | AV18 | 2 | 0.009795, 0.438761 | −1.867744e−4 @ 0.29861 | 7.807344e−4 @ 0.734385 | 0.719432 | 2.09620e−4 |
| CT18NLO | Eq. 21 | **CD-Bonn** | 2 | **0.043538**, 0.391923 | −1.316167e−4 @ 0.27203 | 9.378742e−4 @ 0.729070 | **0.864233** | 4.38973e−4 |
| MSTW | κ = 1 | AV18 | 2 | 0.027946, 0.507853 | −2.083474e−4 @ 0.35706 | 5.646694e−4 @ 0.755642 | 0.520332 | 5.79583e−5 |
| MSTW | κ = 1 | **CD-Bonn** | 2 | 0.064244, 0.465608 | −1.586273e−4 @ 0.33049 | 6.764284e−4 @ 0.745013 | **0.623316** | 2.06925e−4 |
| **MSTW** | **Eq. 21** | AV18 | 2 | 0.027886, 0.495239 | −2.346063e−4 @ 0.35706 | 9.150962e−4 @ 0.771585 | **0.843243** | 2.24896e−4 |
| **MSTW** | **Eq. 21** | **CD-Bonn** | 2 | 0.064129, 0.457018 | −1.769065e−4 @ 0.33049 | **1.085577e−3 @ 0.766271** | **1.000338** | 4.48580e−4 |
| *digitized CDKS* | — | — | *2* | *0.065645, 0.457177* | *−1.768140e−4 @ 0.33236* | *1.085210e−3 @ 0.765663* | *1* | *4.59200e−4* |

The four AV18 rows reproduce §§0–4's and item 4's pinned values to every
printed digit, which is the control that says the CD-Bonn rows were measured
with the same machinery and not a new one.

**The wave-function factor is nearly PDF-independent**, which is the check that
0.843 → 1.000 is not one lucky configuration. AV18 → CD-Bonn on the G3b ratio:

| configuration | factor |
|---|---|
| ToyF2, Eq. 21 | ×1.20150 |
| CT18NLO, Eq. 21 | ×1.20128 |
| MSTW, κ = 1 | ×1.19792 |
| **MSTW, Eq. 21** | **×1.18631** |

A 1.3 % spread across two PDFs and both δ-functions. **§5's residual after A2
was a factor 1.186 unattributed** — and the wave function supplies exactly
1.18631 of it. §5's `item 5` row, which carried the *proxy*'s ×0.88 in the
"opens" column, is superseded: with the real wave function the factor is
×1.186 and it **closes**, which is the sign reversal `phase_A_cdbonn.md` §8
predicted and `phase_D_gate.md:330-331` got backwards.

### 8.3 What got worse, stated rather than buried

* **ToyF2 loses both sign changes.** With CD-Bonn the toy F₂ curve has *no*
  zero in (0, 1] — it is positive from 0.001 up. That is a property of ToyF2's
  shape (both real PDFs put two zeros back, close to the digitized ones), which
  is why checklist item 5 records that row and deliberately does **not** apply
  G3a's counting clause to it. It is also a reason the layer-3 case must stay
  on AV18: switching the always-available row to CD-Bonn would make G3a abort
  in a build with no optional tier.
* **On CT18NLO the picture is mixed.** The second zero moves 0.438761 →
  0.391923, i.e. from missing the digitized 0.457177 by 0.018 to missing it by
  0.065; the dip goes from 5.6 % high to 25.6 % low; the peak position from
  missing by 0.031 to missing by 0.037. Only the peak magnitude and the low-x
  zero improve. **CD-Bonn is not a uniform improvement on an arbitrary PDF.**
* **G3c is still not enforced and still not exact**: 4.48580e−4 against the
  digitized 4.59200e−4 is 0.977 of it (AV18 was 0.490). Much better, but §4's
  point stands — the peak agreeing does not make the integral agree, and here
  the integral agrees too, which is a *new* fact and not a re-statement of the
  peak one.

### 8.4 It is not the sampling grid

The analytic form is sampled on `fdeut.av18`'s own 0.1 fm⁻¹ spacing precisely
so that an AV18 row and a CD-Bonn row differ in the wave function alone.
Refining it (MSTW, Eq. 21, full 300-point scan):

| dk [fm⁻¹] | rows | peak | G3b | zeros |
|---|---|---|---|---|
| 0.10 (the default) | 201 | 1.0855772e−3 | 1.000338 | 0.064129, 0.457018 |
| 0.05 | 401 | 1.0860770e−3 | 1.000799 | 0.064129, 0.457018 |
| 0.02 | 1001 | 1.0861025e−3 | 1.000822 | 0.064129, 0.457018 |

**4.8e−4 relative between the coarsest and the finest, converged by 0.05**, and
the zeros do not move at all. The coarse-grid spline `b1_nuclear.hpp` warns
about is not what produces 1.0003. (The doctest checks this at the peak's own
x rather than over the whole scan, because each finite-|q⃗| point costs ~0.24 s.)

**What this does NOT test** is the parameterisation's fidelity in momentum
space above ≈ 3 fm⁻¹. Machleidt quotes the fit quality only as an *r*-space L²
norm; the ansatz has u ∼ p⁻⁴, w ∼ p⁻⁶ tails which are a property of the ansatz,
not of CD-Bonn, and its ⟨T⟩ = 15.6 MeV has no published CD-Bonn value to check
against (`phase_A_cdbonn.md` §§6, 10). Grid refinement tests the *sampling*, not
the *ansatz*.

### 8.5 What was added, and what it cost

| file | change |
|---|---|
| `include/lipolgen/cluster.hpp` | `CD_BONN_*` (the twenty published numbers, γ, m₀, n, B_d), `struct CdBonnWave`, `cdbonn_wave()`, `cdbonn_fdeut_table()` |
| `src/core/cluster.cpp` | the constraint solve (a Vandermonde/Lagrange closed form, no matrix solve), ψ_S/ψ_D, the closed-form moments, the `FdeutTable` adapter |
| `include/lipolgen/b1_nuclear.hpp` | `enum class DeuteronWaveSource`; `Options::wave` (default `kFdeutFile`), `cdbonn_k_max_fm`, `cdbonn_dk_fm` |
| `src/core/b1_nuclear.cpp` | one branch in the constructor; both paths hand the same `FdeutTable` contract to the same two lines |
| `python/bindings.cpp` | `CdBonnWave`, `cdbonn_wave`, `cdbonn_fdeut_table`, `DeuteronWaveSource`, the three `Options` fields |
| `tests/test_cluster.cpp` | the coefficient gate — needs **no data file**, the parameterisation is analytic |
| `tests/test_b1_nuclear.cpp` | checklist item 5; `GateRow`/`gate_row`/`check_g3a` hoisted out of item 4 so both cases measure with one implementation |
| `python/tests/test_b1_model.py` | P5, four tests |

| suite | before A5 | after A5 |
|---|---|---|
| `./build/lipolgen_tests` | 368 cases / 17206076 assertions / 1 skipped | **370 / 17207963 / 1 skipped, 0 failed** |
| `python -m pytest python/tests -q` | 236 passed | **240 passed** |
| doctest wall clock | 1 m 42 s | 2 m 15 s (item 5 is three threaded 300-point finite-\|q⃗\| scans) |

`b1_nuclear` T9, the 600-row rtol-1e−12 gate against `validation/reference/*.json`,
passes unchanged — the default path is untouched, and a Python test asserts
that setting `wave` explicitly to its default is a bit-for-bit no-op.

### 8.6 Still open after A5

1. **The default is still AV18**, deliberately. The numbers above are a strong
   argument for making `kCdBonn` the default of `DeuteronConvolutionB1` — that
   object exists to reproduce CDKS's figure, and CDKS used CD-Bonn — but doing
   it would move every pinned gate number in layer 3 and item 4 and every
   published row in §§0–7, so it is a decision for a task that says so, not a
   side effect of this one. It is a one-line change plus a re-pin.
2. **The item-5 rescaling proxy still has no committed implementation.** It
   never had one (`grep` finds none); the 2026-09-02 work built it in scratch
   by writing an `fdeut`-format file, and it stays exactly that reachable
   through `Options::fdeut_path`. `phase_A_cdbonn.md` §7 records the two scale
   factors that reproduce it.
3. **The prose that still says the gate fails** is unchanged by A5, including
   `b1_nuclear.hpp`'s own file comment; §7.1's list is untouched and now has
   one more reason to be reconciled.
4. **The e-print vs journal table numbering** (`phase_A_cdbonn.md` §1): the
   coefficients are Table **XX** of `arXiv:nucl-th/0006014`, LaTeX label
   `tab_dwpar`, not the "Tables XVII/XVIII" several planning documents say.
   The published PRC 63 024001 page proofs were not seen. The header cites the
   LaTeX labels, which are unambiguous.
5. **Swapping in one of CDKS's three inputs does not make this their
   calculation.** They used CD-Bonn *and* MSTW *and* their own kernel; the
   kernel here is an independent implementation of their Eqs. (16)/(17)/(21).
