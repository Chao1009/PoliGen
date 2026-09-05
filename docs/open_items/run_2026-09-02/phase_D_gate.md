# Phase D — the A = 2 validation gate, measured

> # ⚠ SUPERSEDED ON G3b — 2026-09-03 (phase A of the next run)
>
> **This document's verdict is out of date. The gate PASSES for the MSTW2008 LO
> nucleon input at CDKS Eq. (21)'s δ-function (`--b1-unpol mstw`), and not for
> the shipped `ToyF2` default, which stays at 0.440 — outside the acceptance
> window and not covered by the lift.** Everything below
> was measured with the library stand-in `ToyF2` (and a CT18NLO proxy for
> checklist item 4, and a *rescaling* proxy for item 5), and on that footing
> G3b was a factor 2.27 low. Both proxies have since been replaced by the real
> inputs:
>
> | row | this document | 2026-09-03 |
> |---|---|---|
> | nucleon PDF (checklist item 4) | CT18NLO stand-in, ratio 0.7193 | **MSTW2008 LO, CDKS's own — ratio 0.843243** |
> | wave function (checklist item 5) | AV18 rescaled to P_D = 4.85 %, ×0.88 (**opens**) | **the real CD-Bonn parameterisation, ×1.18631 (closes)** |
> | G3b verdict | **FAIL** (0.4400 at the ToyF2 default) | **PASS** (0.843243; 1.000338 with CD-Bonn as well) |
> | G3a counting window | [0.02, 1.0] with an unstated scan floor of 0.01 | **(0, 1.0]**, floor **0.001**, both stated in the design |
> | the ⁶Li publication ban | in force | **lifted** (design §5.4 Escalation not triggered) |
>
> **What in this document is still true:** every layer-0/1/2 result, every
> analytic identity, the κ = 1 vs Eq. (21) comparison, the target-mass and R
> rows, and every number labelled `ToyF2` — 0.4400 is still exactly what the
> library stand-in gives. That is a statement about the SHIPPED DEFAULT, not
> about a build option: `ToyF2` is what `--b1-unpol toy` selects and what every
> run that does not ask for another backend uses, in every build, PYTHIA tier
> or not. A build with no PYTHIA tier additionally cannot *reach* the passing
> configuration at all. What is **not** true any more is the verdict, the
> escalation, and the ban — each for the MSTW configuration only.
>
> The current measurements are `docs/open_items/run_2026-09-03/phase_A_numbers.md`
> (§§0–5 the MSTW rerun, §8 CD-Bonn) and `phase_A_cdbonn.md` (the
> coefficients). This file is kept as the record of what was measured on
> 2026-09-02/03, not as a statement of the gate's state.

Design: `docs/open_items/run_2026-09-02/design_D_b1_li6.md` §5.
Code: `include/lipolgen/b1_nuclear.hpp`, `src/core/b1_nuclear.cpp`.
Tests: `tests/test_b1_nuclear.cpp` (T1, four `TEST_CASE`s).
Author: Agent A of design D §8, 2026-09-03; revised the same day after review.

Everything below is **computed by the code in this tree**, never typed.
Reproduce with `./build/lipolgen_tests -tc='b1_nuclear*' -s` (which *pins*
every number except the two that need LHAPDF or a wave function this tree does
not have — checklist items 4 and 5) and with
`python3 validation/b1_li6_table.py`.

**Two x grids appear below and they are not the same one.** The doctest scans
`linspace(0.01, 1.59, 300)`; `validation/b1_li6_table.py` scans the union of
`x_grid(0.01, 1.0, 100)` and `x_grid(0.10, 0.80, 100)` = 200 points, denser
inside the G3b window. Landmarks from the two differ in the 4th digit for that
reason alone (second zero 0.37736 against 0.37737, G3b ratio 0.43999 against
0.43997). **Every table below says which grid made it.**

**One convention for the G3b number.** The ratio *R* = max|x·b₁| / the
digitized peak and "a factor 1/*R* low" are the same statement; both are
quoted, always together, never one alone.

> **REVISION OF 2026-09-03 (review).** `DeuteronConvolutionB1::Options::
> finite_q_delta` now defaults to **true**, i.e. CDKS Eq. (21)'s exact
> δ-function y = (E − p_z κ)/M_N with κ = |q⃗|/ν, and every gate number below
> is quoted there. Eq. (17)'s (E − p_z)/M_N is the "≃" of their Eq. (18) — the
> κ → 1 limit — so quoting the *gate*, whose whole job is to reproduce CDKS's
> own figure, at κ = 1 overstated the deficit by the design's own checklist
> item 0: 3.68 rather than 2.27. The verdict does not change (2.27 is still
> outside G3b's factor-2 window). `Li6ConvolutionOptions::finite_q_delta`
> keeps its `false` default (design A10): in ⁶Li the same switch is worth
> −1 % to +7 %.

---

## Verdict

| clause | result |
|---|---|
| **G0a–G0g** (layer 0, analytic) | **PASS** |
| **G1a–G1f** (layer 1, unpolarized convolution) | **PASS** |
| **Layer 2** (α–d sign gate) | **PASS** |
| **G3a** (shape, HARD) | **PASS** at the default (κ = √(1+γ²), ToyF2, `r1998`) — but see the margin note below |
| **G3b** (magnitude, within a factor 2) | **FAIL** at the default: ratio **0.4400**, i.e. a factor **2.27 LOW** (0.2717 / factor 3.68 at Eq. (17)'s κ = 1) |
| **G3c** (Close–Kumano, reported) | recorded below |

**The gate as a whole therefore FAILS**, on G3b only, and design §5.4's
*Escalation* clause applies: the backend stays in the tree behind its flag with
a `WARNING: A = 2 gate not passed` in the header, and **no ⁶Li number may be
published** while that stands. Item 10 stays open.

> **SUPERSEDED 2026-09-03.** The paragraph above is the verdict *at the ToyF2
> default*, and it is still what that configuration gives. It is no longer the
> gate's verdict: with checklist item 4 worked (MSTW2008 LO, CDKS's own PDF)
> G3b is **0.843243** — inside the factor-2 window — so Escalation is not
> triggered, the `WARNING` is gone from the header and `--help`, and **the ⁶Li
> publication ban is lifted** — for that configuration (`--b1-unpol mstw` at CDKS
> Eq. (21)); on the shipped `toy` default the same clause gives **0.440**,
> outside [0.5, 2], so the default's numbers are **not covered by the lift**.
> The mandatory ±100 % band is not lifted either.

**G3a passes, but read how.** The low-x zero clears its counting window's
lower edge (x = 0.02) by **0.0021**, and with a realistic PDF (CT18NLO,
checklist item 4) it drops below the scan floor altogether, so the "exactly two
sign changes in [0.02, 1.0]" clause *fails* there. **The low-x zero is not a
robust discriminator.** The second zero (0.3774 against the digitized 0.4572,
inside its ±0.10 window by 0.020) and the peak position (0.755 against 0.766)
are.

**What closes G3b.** The remaining factor 2.27 is attributed to one identified
item: the **nucleon PDF** (item 4). With `LhapdfSF("CT18NLO")` in place of
`ToyF2` the ratio is **0.7193** (a factor **1.39** low), *inside* G3b's
factor-2 window. It is not this phase's default to change — the core must not
link LHAPDF, and CDKS used MSTW2008 LO, which is not installed — so the finding
is recorded, not absorbed.

> **RESOLVED 2026-09-03, and the last clause of that paragraph was wrong.**
> MSTW2008 LO *is* on disk — the CENTRAL member ships with PYTHIA 8 as
> `pdfdata/mstw2008lo.00.dat`, read by `Pythia8::MSTWpdf`; only LHAPDF's set
> store lacks it. `MstwSF` (`include/lipolgen/mstw_sf.hpp`) reads it with the
> same charge weights as `LhapdfSF`, and with it the ratio is **0.843243**,
> not 0.7193. The prediction in this paragraph — that the nucleon PDF closes
> G3b — held; the stand-in simply understated it by ×1.17.

---

## Layer 0 — analytic identities (no data, no figure)

| id | statement | measured | verdict |
|---|---|---|---|
| G0a | ⟨P_zz⟩_{L=2} = +1/10 (H = ±1), −1/5 (H = 0); δ_T ratio D-wave/S-wave = (−3/10)/(−3) | exact rational arithmetic in the test | **exact** |
| G0b | ⟨S_z⟩_{L=2} = −1/2 (H = +1) | exact rationals; matches `ALPHA_D_VECTOR_POLARIZATION` = 1 − (3/2)P_D | **exact** |
| G0c | Eq. (21) SD coefficient −3/(4√2π) from 2·Y₀₀Y₂₀·(C⁰₀ − (C⁺¹₀+C⁻¹₀)/2) | −0.1688094 vs −0.1688094 | 1e-12 |
| G0d | Eq. (21) DD coefficient +3/(16π) from (3/10)\|Y₂₀\|² + (3/10)\|Y₂₁\|² − (3/5)\|Y₂₂\|² | identity in c², checked at five c | 1e-12 |
| G0e | ∫δ_T f dy = 0 and ∫(dy/y)δ_T f = 0 | see table below | see below |
| G0f | P_zz-weighted S–D CG sum = +0.9486833 = −(Eq. 21's), so the dropped b₁ᵈ-sector SD coefficient is exactly 1/3 of Eq. (21)'s | 1/3 exactly | 1e-12 |
| G0g | `k_range` against a brute-force \|c*\| ≤ 1 scan, κ = 1 and 1.38, all three mass assignments | both endpoints to the scan's own step (<= 8e-5), \|c*\| = 1 at both endpoints to 1e-9 | **PASS** |

**G0e in detail.** The identity is exact for *any* wave function and *any* k
cutoff, but only over the **whole** y support, which reaches **negative** y —
E = M − ε − k²/(2M_recoil) turns negative above k = √(2 M_s M_r), and
y = (E − kκc)/M_s with it. The design's grid floor y_min = 10⁻⁴ truncates that
tail. Measured (residual as a fraction of max\|δ_T f\|):

| density set | y_min = 10⁻⁴ (the default) | y_min = −3 |
|---|---|---|
| deuteron, ∫δ_T f dy | −1.129e−4 | **−4.2e−7** |
| deuteron, ∫(dy/y)δ_T f | +4.53e−4 | **+7.5e−8** |
| ⁶Li struck d, ∫δ_T f dy | +4.9e−8 | +7.7e−7 |
| ⁶Li struck α, ∫δ_T f dy | −2.6e−10 | −6.3e−9 |

The deuteron residual on the default grid does **not** improve with grid
refinement (1.129e−4 → 1.129e−4 → 1.129e−4 at ×1/×2/×4) — it is a *support*
truncation, not a discretization error. The α–d sets never see it, because
their k table stops at 5 fm⁻¹ where E is still positive. `y_min` is an
`Options` field (default unchanged) precisely so this can be shown.

---

## Layer 1 — the unpolarized convolution

AV18 u(k), w(k) from `data/vmc/deuteron/fdeut.av18` with the file's own
ε_d = 2.224574 MeV, before renormalisation.

| id | quantity | design | measured | verdict |
|---|---|---|---|---|
| G1a | ∫k²[u²+w²]dk, the file's own 0.1 fm⁻¹ grid | 0.999976 | **0.9999764** | PASS (1e-4) |
| G1a′ | the same, **spline**-refined on 10⁵ points | — | **1.0004598** | recorded |
| G1a″ | the same, **linear**-refined | (design §3: 1.0218) | **1.0217771** | the 2.13 % bias the spline removes |
| G1b | P_D of the file | 0.0575999 | **0.0575999** | PASS (1e-5) |
| G1c | ∫f dy = 1 − ε/M_N − ⟨p²⟩/(2M_N²) | 0.98707 | **0.9876619**; the moment identity from the same table gives 0.987066 | PASS (1e-3), 6.0e-4 apart |
| G1d | mean_y = ⟨y²⟩_φ/⟨y⟩_φ | 0.99546 | **0.9952583**; identity 0.995454 | PASS (1e-3) |
| G1e | y_max at κ = 1 | 1.4976 | **1.497629** | PASS |
| G1e | y_max at x = 0.8, Q² = 2.5 | 1.9502 | **1.94840** | PASS (9.2e-4, just inside 1e-3) |
| G1f | F₁ᴰ/F₁ᴺ per nucleon at x = 0.1/0.3/0.5/0.8 | 0.997/0.990/0.986/1.108 | **0.99754 / 0.99162 / 0.98429 / 1.07082** | first three PASS at 1 %; **x = 0.8 is 3.4 % off** and is pinned to the measured value |

The Fermi-motion/EMC shape is there: a dip of 1.6 % near x = 0.5 and a rise
above x ≈ 0.65. The design's G1f row does not say which R or which target-mass
convention its x = 0.8 entry was computed with; the three low-x points agree at
the stated 1 %, so the kernel is not what differs.

**Note on G1c and the interpolation (checklist item 6).** ∫f dy is
*necessarily* < 1 with CDKS Eq. (12) — the p_z term averages to zero. The
linear-interpolation artefact the design warns about is reproduced exactly:
1.0217771/1.0004598 = 1.021308, and the baryon-number renormalisation would hide it.
The wave functions are splined (`ClusterPartialWave::operator()`), and G1c is
confirmed *before* anything downstream is believed.

---

## Layer 2 — the α–d sign gate

| statement | measured | verdict |
|---|---|---|
| Q of the α–d relative motion < 0 | **−0.3333 fm²** | **PASS** |
| the same code on the deuteron pair reproduces `fdeut.av18`'s own header `qm` = 0.269673 fm² | **0.269362 fm²** (0.12 %) | **PASS** — this validates the machinery, not just its sign |
| φ₀φ₂ (CDKS convention) **> 0** at k = 0.2 fm⁻¹ (below the S node) | +5.87e+1 | **PASS** |
| φ₀φ₂ **< 0** at k = 1.0 fm⁻¹ (between the nodes) | −9.60e+0 | **PASS** |
| φ₀φ₂ **> 0** at k = 3.0 fm⁻¹ (above the D node) | +9.52e−4 | **PASS** |
| the **deuteron's** φ₀φ₂ < 0 at both 0.2 and 1.0 fm⁻¹, same code path | −1.22e+2, −8.41e+0 | **PASS** |

So the α–d S–D relative sign is opposite to the deuteron's below the S node
(0.678 fm⁻¹, 66 % of the density) and deuteron-like between the nodes, exactly
as design §2.2's node table states, and the quadrupole built from the same two
waves comes out with the sign of the measured Q(⁶Li).

---

## Layer 3 — CDKS Fig. 4

`b1_landmarks_of_table()` on the raw digitized column, **computed**:

| landmark | this code |
|---|---|
| sign change, falling | **0.065645** |
| minimum | **−1.76814e−4** at x = **0.33236** |
| sign change, rising | **0.457177** |
| maximum | **+1.08521e−3** at x = **0.76566** |
| ∫b₁ dx over the digitized range | **+4.5920e−4** |

which reproduces every entry of design §5.4's table, and equals
`close_kumano_integral(true)` to 1e-9.

### The kernel, DEFAULT options (ToyF2, R = `r1998`, κ = √(1+γ²))

**Grid: the doctest's `linspace(0.01, 1.59, 300)`** (`tests/test_b1_nuclear.cpp`,
"T1 gate layer 3"), which pins every number in this table.

| landmark | this kernel | the κ = 1 form | digitized | Δ (default) | G3a window |
|---|---|---|---|---|---|
| first sign change (falling) | **0.022113** | 0.022137 | 0.065645 | 0.0435 | ±0.08 → **PASS** |
| second sign change (rising) | **0.377357** | 0.392321 | 0.457177 | 0.0798 | ±0.10 → **PASS** (by 0.020) |
| peak position in [0.5, 1.0] | **0.75508** | 0.73923 | 0.76566 | 0.0106 | ±0.10 → **PASS** |
| number of sign changes in [0.02, 1.0] | **2** | 2 | 2 | — | **PASS** |
| minimum | −3.30912e−5 at 0.23722 | −3.12115e−5 at 0.23722 | −1.76814e−4 at 0.33236 | — | not gated |
| maximum | +4.77477e−4 at 0.75508 | +2.94842e−4 at 0.73923 | +1.08521e−3 at 0.76566 | — | G3b below |
| ∫b₁ dx, [0.01, 1.59] | **+2.15370e−4** | +1.07915e−4 | +4.5920e−4 | — | G3c |

> **G3a's two fragile clauses, flagged not hidden.**
> (a) The first sign change sits at 0.022113, i.e. **0.0021 above the
> [0.02, 1.0] counting window's lower edge**. A 10 % move of the low-x tail
> turns "exactly two sign changes in [0.02, 1.0]" into one and G3a fails on a
> *counting* clause rather than on physics — and that is **not hypothetical**:
> with CT18NLO (item 4) the zero drops below the scan floor and the clause does
> fail. The doctest asserts `z[0] > 0.02` and pins the margin explicitly, so
> the failure, if it comes, is legible.
> (b) The second zero clears its ±0.10 window by **0.020**.
> **The low-x zero is not a robust discriminator; the second zero and the peak
> position are.**

**G3b.** max|x·b₁| over [0.10, 0.80] = **4.77477e−4** against the raw column's
1.08521e−3: ratio **0.4400**, a factor **2.27 LOW**. Outside the factor-2
window ⟹ **FAIL**. (At Eq. (17)'s κ = 1: 2.94842e−4, ratio **0.2717**, a factor
**3.68** — the number this document quoted before the review of 2026-09-03.)

**G3c.** ∫b₁ dx, this kernel **+2.15370e−4** (κ = 1: +1.07915e−4);
`close_kumano_integral(true)` (digitized CDKS) **+4.59200e−4**;
`close_kumano_integral(false)` (Miller) **+5.91474e−3**. None is zero; none is
enforced.

---

## The §5.4 checklist, item by item

### 0. The δ-function — κ = 1 (Eq. 17) or κ = √(1+γ²) (Eq. 21)?

**Run both ways, as design §5.4 requires.** Ratio κ/(κ = 1) of x·b₁, doctest
grid:

| x | κ = 1 (Eq. 17) | κ = √(1+γ²) (Eq. 21) | ratio | design's expectation |
|---|---|---|---|---|
| 0.10 | −1.99274e−5 | −2.03395e−5 | 1.021 | 1.02 |
| 0.20 | −3.06060e−5 | −3.22675e−5 | 1.054 | 1.05 |
| 0.30 | −2.78408e−5 | −2.85688e−5 | 1.026 | 0.94 |
| 0.50 | +8.44989e−5 | +1.32595e−4 | **1.569** | 1.54 |
| 0.70 | +2.85771e−4 | +4.51252e−4 | **1.579** | 1.59 |
| 0.80 | +2.74790e−4 | +4.63725e−4 | **1.688** | 1.67 |
| ∫b₁ dx (80-pt grid) | +1.10438e−4 | +2.17895e−4 | **1.973** | 1.94 |

Landmarks with κ (doctest grid): zeros 0.022113 (falling) and 0.377357
(rising), peak +4.77477e−4 at 0.75508. **G3a still passes with κ**
(|0.377357 − 0.457177| = 0.0798 < 0.10), and G3b improves from **0.2717** to
**0.4400** — a factor 3.68 → 2.27, still outside.

*The design's revision-2 note asked the implementer to measure, not assume, the
zero shift. Measured: κ moves the second zero from 0.392321 **down** to
0.377357, i.e. slightly **away** from the digitized 0.457177, which is the sign
this session's design re-run predicted and the opposite of the adversarial
review's. The magnitude effect is confirmed.*

**DEFAULT CHANGED, 2026-09-03 (review): κ = √(1+γ²) on this object.** Eq. (21)
is CDKS's exact δ-function — y = M p·q/(M_N P·q) = (E − p_z κ)/M_N, their
Eq. (18) — and Eq. (17)'s (E − p_z)/M_N is the "≃" of it. A gate whose whole
job is to reproduce *their* figure is quoted at *their* equation; design §5.4's
own item 0 says a gate-driven change here "is a finding, not a fudge". Two
further facts point the same way: at κ = 1 the kernel's support ends at
y_max = 1.4976, so x·b₁ ≡ 0 above x = 1.4976 while the digitized theory-1
column is still +4.04e−6 at x = 1.590 and CDKS's text stresses "even extremely
large x (x > 1)"; and the κ = 1 column overstated the deficit by 62 %
(3.68 against 2.27).
`Li6ConvolutionOptions::finite_q_delta` is a different object and keeps its
`false` default — in ⁶Li the switch is worth −1 % to +7 % (design A10).

### 1. Target mass — is CDKS Eq. (22) in?

Yes, and at the default δ-function it is worth **×1.48729** at x = 0.8
(4.63725e−4 with, 3.11791e−4 without); at κ = 1 it is **×1.52049**.
`NuclearF2::f1a` appears nowhere in `b1_nuclear.cpp`; T13 pins the ratio of the
kernel's own F₁ᵈ to `f1a`/2 at exactly 1 + γ².

### 2. R — `r1998`, not `r_sigma_lt`

`r1998` (CDKS's SLAC world fit) is the **default of the gate**. Switching to
the programme's toy R moves the second zero from **0.392321 to 0.365215**, i.e.
*away* from the digitized 0.457177, and the peak from 2.94842e−4 at 0.73923 to
2.99437e−4 at 0.75508 (+1.6 %). Worth a few % as the design says, and it is the
choice that makes the *shape* better.

**Measured at κ = 1, deliberately, and that is the whole of the difference from
the table above.** R is a property of F₁ and is independent of the
δ-function, while a 300-point landmark scan at the default κ costs 73 s (the
densities are rebuilt at every x); the r1998 column it is compared against is
the κ = 1 one, so the two are like for like. The doctest pins both.

**The ⁶Li backend does NOT default to `r1998`.** `Li6ConvolutionOptions::r_func`
null means `r_sigma_lt`, the kernel's own toy R, for consistency with
`InclusiveKernel`'s F₁ — while this gate ran with `r1998`, CDKS's own. The
difference is the +1.6 % peak and the 0.392 → 0.365 second zero measured just
above, i.e. *not* cosmetic but far inside the 100 % band. It is stated in
`docs/USAGE.md` §2a and `OPEN_ITEMS_SOLUTIONS.md` §10; whether the ⁶Li default
should become `r1998` is a physics decision, not a bug, and it is an **open
follow-up, not decided in the 2026-09-03 close-out** (`OPEN_ITEMS_SOLUTIONS.md`
§10 item 4).

### 3. Per nucleon vs per deuteron — RESOLVED, not re-litigated

The gate compares against the **raw column** and the ⁶Li backend's default
b₁ᵈ is `cdks_b1_raw_per_nucleon()`, which is `tables::kB1CdksQ2p5()` without
`B1_PER_DEUTERON_TO_PER_NUCLEON`. That constant is **not** changed
(`sf.hpp`/`constants.hpp` are on this phase's do-not-touch list); Q1b stays a
separate close-out item.

### 4. Nucleon PDF — ToyF2 vs a real set: **THE DOMINANT ITEM**

Rerun with `LhapdfSF("CT18NLO")` from outside the core (the core must not link
LHAPDF).

**Grid: `validation/b1_li6_table.py`'s own** — the union of
`x_grid(0.01, 1.0, 100)` and `x_grid(0.10, 0.80, 100)`, 200 points, denser
inside the G3b window — **not** the doctest's `linspace(0.01, 1.59, 300)` of
the Layer 3 table above. That, and nothing else, is why the two disagree in the
4th digit (0.37737 here against 0.37736 there). Reproduce with
`python3 validation/b1_li6_table.py`:

| configuration | second zero | minimum | peak | G3b ratio |
|---|---|---|---|---|
| ToyF2, κ (**the default**) | 0.37737 | −3.3089e−5 @ 0.2343 | 4.7751e−4 @ 0.7576 | **0.4400** (×2.27) |
| ToyF2, κ = 1 | 0.39230 | −3.1211e−5 @ 0.2343 | 2.9484e−4 @ 0.7400 | **0.2717** (×3.68) |
| CT18NLO, κ | 0.43877 | −1.8679e−4 @ 0.3000 | **7.8062e−4 @ 0.7364** | **0.7193** (×1.39) |
| CT18NLO, κ = 1 | **0.44920** | **−1.7025e−4 @ 0.2980** | 4.9178e−4 @ 0.7152 | **0.4532** (×2.21) |
| *digitized CDKS* | *0.45718* | *−1.76814e−4 @ 0.3324* | *1.08521e−3 @ 0.7657* | *1* |

A real PDF set fixes **both** things ToyF2 gets wrong: the second zero lands at
0.449 against the digitized 0.457, and the mid-x dip comes out −1.70e−4 against
−1.77e−4 — a 4 % agreement where ToyF2 was a factor 5.7 low. On the default
δ-function the peak ratio is **0.7193**, a factor 1.39 low: *inside* G3b's
window. **This is now the only identified item still open at the peak.**

*Caveat, stated: CDKS used MSTW2008 LO, not CT18NLO. CT18 is a stand-in of the
right kind (a real global fit with a realistic large-x valence exponent), not a
reproduction of their input.*

**One clause regresses.** With CT18 the first sign change moves below x = 0.01
(the scan floor), so "exactly two sign changes in [0.02, 1.0]" becomes one and
G3a's counting clause fails. That is the same fragility flagged above, seen
from the other side: the low-x zero position is not a robust discriminator.

### 5. Wave function — AV18 (5.76 %) vs CD-Bonn (4.85 %)

Crude proxy: AV18's w(k) rescaled so that P_D = 4.850 %, S renormalised,
measured **at κ = 1** (before the default flip; the ratio, not the third digit,
is the finding, and no CD-Bonn table exists in this tree to do better). Peak
**2.5855e−4 at 0.7392** against AV18's 2.9437e−4 — i.e. **12 % LOWER**, moving
*away* from the digitized curve, and the second zero moves to 0.40952 (towards
it). The residual factor is therefore **not** explained by the wave function;
if anything CDKS's softer D state should make their curve smaller, not larger.
*This proxy rescales the D-state normalisation only and does not reproduce
CD-Bonn's softer high-k tail, which CDKS say "plays an important role"; a real
CD-Bonn table would be needed to settle it.*

### 6. Interpolation — confirmed before anything else

`∫k²(u²+w²)dk` = 0.9999764 on the file's own grid, **1.0004598** splined on a
fine grid, **1.0217771** linearly interpolated on the same fine grid. G1c holds on
the spline (0.98766 against the moment identity's 0.98707). The 2.13 % linear
bias, and the fact that the renormalisation hides it, are both reproduced.

---

## Residual budget

Starting point: the κ = 1 column, ratio 0.2717 (a factor 3.68 low).

| piece | measured effect on the peak | direction |
|---|---|---|
| finite-\|q⃗\| δ-function (κ = 1 → √(1+γ²)) | ×1.620 | closes — **now IN** (the default) |
| nucleon PDF (ToyF2 → CT18NLO) | ×1.669 | closes — **still open** |
| target mass (already in) | ×1.49 at κ (1.52 at κ = 1) | in |
| R (`r_sigma_lt` → `r1998`, already in) | ×0.98 | in |
| wave function (AV18 → P_D = 4.85 %) | ×0.88 | **opens** — **SUPERSEDED**: that is the *rescaling proxy*. The real CD-Bonn is **×1.18631** and **closes** (`phase_A_numbers.md` §8) |
| figure digitization | not quantified | — |
| **the default today** | **ratio 0.4400 (a factor 2.27 low)** | **outside G3b** |
| **with the PDF item as well, measured** | **ratio 0.7193 (a factor 1.39 low)** | **inside G3b** |

The two are very nearly independent: taken separately they predict
0.2717 × 1.669 × 1.620 = **0.735** against the measured **0.7193**, a 2 %
overlap, so there is no third unexplained factor at the peak. **What remains
unexplained is the low-x tail**, where even CT18 + κ gives
x·b₁(0.10) = −3.8e−5 against the digitized −1.73e−5.

## What a follow-up must do

1. ~~Decide item 0's default~~ — **done, 2026-09-03**: the gate defaults to
   CDKS Eq. (21)'s δ-function. The ⁶Li backend's own default stays Eq. (17)
   (design A10), where the switch is worth 1 %.
2. Get MSTW2008 LO, or ask S. Kumano's group for the tabulated curves (design
   Q2), and rerun item 4 with CDKS's own PDF. **This is now the only
   identified open item at the peak.**
3. Get a real CD-Bonn u, w rather than the rescaling proxy of item 5.
4. Decide whether `Li6ConvolutionOptions` should default to `r1998` like the
   gate rather than to `r_sigma_lt` (item 2).
5. Only then re-open G3b. Until it passes, **no ⁶Li number ships.**

> **DONE — 2026-09-03.** (2) `MstwSF` reads the MSTW2008 LO grid PYTHIA ships;
> G3b 0.7193 → **0.843243**. (3) The real CD-Bonn Appendix-D parameterisation
> is `cdbonn_wave()` in `cluster.hpp`, reachable as
> `DeuteronConvolutionB1::Options::wave = kCdBonn`; it is **×1.18631** on the
> peak, i.e. it *closes* rather than opens, which is the sign this document's
> rescaling proxy got backwards. (4) Decided: `Li6ConvolutionOptions` keeps
> `r_sigma_lt`, because the tensor weight is a ratio whose denominator carries
> `InclusiveKernel`'s R — recorded with its measured cost in
> `docs/OPEN_ITEMS_SOLUTIONS.md` §10. (5) G3b re-opened and **passes at
> `--b1-unpol mstw`, Eq. (21)** — not at the `ToyF2` default, which stays at
> 0.440; the ban is lifted for the configuration that was measured. The
> residual budget above is superseded by
> `phase_A_numbers.md` §5 and §8.2.
