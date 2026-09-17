# Run 2026-09-03 — one page

What the open-items run measured, phase by phase, in the numbers a reader would
cite — **each with the window it was measured in** — and, at the end, what it
did **not** establish. Every number here is quoted from the phase record named
beside it; `STATUS.md` is the board and the decision registry, `PLAN.md` the
brief, `AUTHOR_DECISIONS.md` the decisions stated in full — **25 when this
page was written, 27 since 2026-09-15**, when rows 26 and 27 were added for
the tagged S–D interference phase and the sibling repository's published
numbers (see item 9).

**Start → end.** Baseline `a94fd6e`: 363 doctest cases / 17 203 863 assertions /
1 skipped, 211 pytest. End of phase F, re-measured 2026-09-05: **401 doctest
cases / 17 240 286 assertions / 1 skipped / 0 failed**, **926 pytest passed /
112 skipped**, strict citation gate **1206 refs / 96 ranges / 7 external / 0
broken / 6 allow-listed** (re-measured 2026-09-16 from `git archive bd775bc`;
**1205 / 97** was published on the day, off the pre-commit working tree) (+ 19/115 on `SPIN32_FINITE_GAMMA.md`, 8 external on
`PYTHIA_BRIDGE.md`), **SPDX 101/101**. All nine `validation/reference/*.json`
are byte-identical to `a94fd6e`: **no shipped default moved in six phases**, and
the rtol-1e−12 gates did not move.

> **Two of the nine have moved since — 2026-09-06, and deliberately** (item 6b
> below; registry row 26 / `AUTHOR_DECISIONS.md` §B27). Re-measured
> 2026-09-15 by sha256 against `a94fd6e`: `b1_default_li6`, `beams`,
> `bookkeeping`, `coherent`, `spectator`, `spin` and `xsec` are **still
> identical**; `tagged.json` differs in its two spin-1 `model` blocks and
> `_manifest.json` in the one provenance line that names it. So the sentence
> above holds **as the record of this run** (`a94fd6e..dffe94e`) and, read as
> a statement about the tree today, must be read as *byte-identical except
> those two files, by row 26*. The follow-on run's own start→end line is in
> `../run_2026-09-06/SUMMARY.md`.

## A — the b₁(⁶Li) publication ban (`ac22331`, `phase_A_*.md`)

* **G3b = 0.843243** — max\|x·b₁\| over x ∈ [0.10, 0.80] against the digitized
  CDKS Fig. 4, with **MSTW2008 LO** at **CDKS Eq. (21)**'s δ-function. Inside
  the factor-2 window; the lower edge is cleared by 1.69.
* **0.440 on the shipped `ToyF2` default** — the same row, *outside* [0.5, 2].
  The lift is about a **configuration**, not a build: quote `--b1-unpol mstw`.
* **0.520** at Eq. (17)'s κ = 1; **1.000338** with CD-Bonn *and* MSTW together —
  a residual below the reference figure's own digitization error, not
  three-digit agreement with CDKS.
* The ban is lifted across the **17 files** that carry it (`git grep -il
  'publication ban'`, the same 17 at `ac22331` and at HEAD). The **±100 % band
  stays**: the gate is A = 2 and says nothing about the α–d step.

## B — the radiative tail (`d3ac125`, `phase_B_numbers.md`)

* POLRAD's "only the t-peak leads" holds **event-weighted (+0.61 %)** in the
  Q² ≥ 20 GeV², y ≤ 0.9 window (5182 of 200 000 events, seed 1234, **at
  `--pz 0`** — P_z named 2026-09-15 and both fills re-measured: the CLI default
  P_z = 0.7 gives 5194 events and +0.62 %; the 5182 here is exact) and **fails
  per cell**: **331 of that window's 1356 accepted cells (24.4 %)** differ by
  > 1 %, worst **×6444** at x = 0.7943, y = 0.0088, **318 of them at y < 0.1**.
* The ⁶Li **C0 shape is a band in both sign and magnitude**: σ^el_T/σ^el_U at
  x = 0.10, Q² = 5 is **+9.357e−04 on the `ho` edge** and **−2.442e−04 on
  `vmc-ft`** (`rc.hpp:1150`/`:1155`; six times §B3.1's `(1/6)σ^el_T/σ^el_U`
  column, +1.559536e−04 / −4.070157e−05, which is the form the phase record
  tabulates) — a factor **68** resp. **≈ 260** below POLRAD's own deuteron
  elastic-tail value +0.064 — so the published **sign** is withdrawn as a band
  edge, and so is the magnitude.
* The polarised quasi-elastic tail gets a **price tag, not a calculation**:
  Δ(ΔA_zz) = **0.086 % / 0.0003 % / 0.25 %** of the band at x = 0.01 / 0.10 /
  0.30, from a **borrowed** magnitude.
* A = 2 → A = 6 transfer of δ(x): A_zz band half-widths **1.358472e−04 /
  1.518818e−04 / 1.921170e−04** at f = 0 / 0.5 / 1, x = 0.01, Q² = 5.
* The low-x anchor's source has **0 INSPIRE citations** (record 647050).

## C — coherent ⁶Li and the VMC inputs (`adec442`, `phase_C_numbers.md`)

* `eps_b0` = −0.08 implies **Q_charge(⁶Li) = −0.9345 fm², 11.42×** the measured
  −0.0818; the honest ⁶Li band at `slope_b` = 50 is **−(0.0070 … 0.0527)**. The
  default is **kept** (a reference gate) with that cost recorded.
* One deuteron per run: the α-tag channel had held the VMC α–d overlap and a
  Hulthén embedded deuteron together — **+2.069 %** on every polarized tagged-α
  observable (ratio **1.020687209533**, exact at every (x, Q²)). Both now follow
  `--cluster-wave`; the opt-in path is **−2.027 %**, the default **bit for bit**.
  *(**Re-measured 2026-09-06** after the tagged sector's S–D interference sign
  was fixed: **+2.069 % survives** — the ratio of the two channels'
  `vector_dilution` is **1.020687624664**, against **1.020687500612** before,
  because both dilutions are angle-integrated and the flipped S–D cross term
  integrates to zero. The **−2.027 %** likewise: −2.02683 % against
  −2.02682 %.)*
* AV18 deuteron **P_D = 0.057600**; a `vmc` whole-nucleus reading is
  **0.887076** (measured end to end at 0.886169 ± 0.0011 on 400 k events)
  against the ab-initio **0.848** and the shipped **0.811228**.
* **O5 is MARGINAL, and it is a band**: the ⁶Li tensor a₂ in coherent J/ψ at one
  EIC year is **S = 2.63 σ** at the low edge with the top a **span 2.84 … 3.29
  σ**, **3 σ at 8.3 … 13.0 fb⁻¹/u** — quoted at 10 × 99.5 on the
  background-immune ⟨P_zz²⟩ = 0.81, from a **closed-form map, not a dipole
  amplitude**. The ×13.5 in rate is Q² < 0.1 (**85 %** of the coherent J/ψ rate)
  plus the μ⁺μ⁻ channel.
* **O3** is bounded, not answered: a committed correlated input (`he4.dd`)
  bounds 4-body correlation on a comparable position-space moment at **+18.2 %
  (variance) / +8.7 % (rms)** against the uncorrelated-product null.
* The citation gate **could not fail** before this run: it accepted any match
  within ±2 lines. Strict now — **160 references re-anchored across 27 files**,
  and today's gate prints **239 broken** on (`31ed181` doc, `31ed181` sources).

## D — the sweep (`de1a040`, `phase_D_*.md`)

* **FSI**: the two Glauber variants are not one — **99.50 %** of events differ
  by more than **1 %** (|w_nucleon/w_cluster − 1| > 0.01), ratio to **68.5**, on
  ONE stream (`tagged-6Li-alpha`, 20 000 events, seed 1234, `tensor-thirds`,
  **σ_XN = 40 mb**, one end of the mandatory 20–40 mb band). Σw/Σw_off
  **0.520954** vs **0.581605**; the model's own grid-integrated `survival()`
  **0.520239** vs **0.582899** — not the same quantity, and four sites had
  called the first the second.
* **PomSet**: over 15 sets × 20 000 events the T0 columns are **bit-identical
  across the fourteen sets that still run**, so that band is identically zero;
  the systematic is hadronic — ⟨n_charged⟩ **−2.6 % / +9.9 %** over the twelve
  DPDF fits, kaons **×2.8**. Set 11 is refused (100 % e_q² fallback).
* **The coherent \|t\| ceiling stays 0.2** and its *reason* changed: the anchor
  range (\|t\| ≤ 0.30, knob-independent) is primary, positivity contingent —
  edge **0.245** at the shipped `eps_b0`, **2.80** at the measured quadrupole.
  On 200 000 events (seed 99) ⟨\|t\|⟩ = **0.019963**, max **0.197575**, and the
  truncation redistributes **4.5e−5** of the rate.
* **Structure-function backends** reach every kernel: cost of the toy on ⁶Li at
  config 1 is σ **×0.7985** (ct18nlo) / **×0.7934** (mstw), A_zz ×1.2524 /
  ×1.2604 — and the toy g₁ⁿ has the **wrong sign over ≈ 0.25 < x < 0.6**.
* **⁷Li**: the rank-2 sector is **exactly 0** and now says so everywhere. The
  α–t convolution is **2.99 ± 0.02 ×** ⁶Li's orbital term *in the research note
  only*; its sign flips with the unpolarised backend, so it is not shipped.
* **Knob provenance is one mechanism**: 66 rows on a default run, driving meta
  and banner, with a matrix pytest of **12 specs × 71 variants = 547 cells**.

## E — release hygiene (`7f68339`, `phase_E_numbers.md`)

* The strict gate covers **three** documents; `SPDX-License-Identifier` on
  **101 files**, gated; `CITATION.cff` and `AUTHORS` exist with a **placeholder**
  where the name goes.
* `ci.yml` (4 jobs) + `build_deps.sh` exercised locally end to end: deps
  **337.28 s** at −j8, job 2 reproducing this machine's **400 cases /
  17 240 262 assertions bit for bit at the commit it was measured against**
  (`7f68339`; phase F then added one case, so the tree now reads 401 /
  17 240 286), cache **431 MiB on disk / 113.5 MiB compressed**. **Never run
  by GitHub.**
* `auditwheel repair`: **1 779 767 B → 7 232 971 B (×4.06)**, tag
  **`manylinux_2_35_x86_64`**. Nothing was uploaded anywhere.

## F — close-out (commit *pending*)

48 review findings applied under five lenses; two touch behaviour — the route
knobs' labels keep the **envelope name** on a channel that has a route (the
classifier *is* consulted there; only the sensitivity is absent), and `--lumi X`
alone now reaches luminosity mode. **Twenty-two** data-gated test cases moved from an
in-case `MESSAGE` (tallied PASSED) to registration-time `doctest::skip`
(tallied SKIPPED) — 8 in `test_tagged.cpp`, 5 in `test_cluster.cpp`, 4 in
`test_mstw_sf.cpp`, 3 in `test_coherent.cpp`, 1 in `test_pipeline.cpp` and
the new `item 5v` — so a case that is nothing but the MSTW grid or `data/vmc`
no longer reads as a pass that measured nothing. **Three in-case guards
survive on purpose**, each inside a case whose shipped-default rows must still
run: T16's MSTW reach row (`test_b1_nuclear.cpp:2034`), T17's verdict row
(`:2147`) and T1r's companion report (`:892`) still print a SKIPPED message
inside a case doctest tallies as passed — measured with the grid absent,
`-tc="*T16*,*T17*"` reports `7 passed | 0 failed`. That is decision row 25,
open (§B26). The close-out then rewrote the four reader-facing documents and
added this one.
It also caught two more instances of the run's recurring defect — a claim
stated more strongly than it was measured: the `RouteReach` header stated a
route-label reach **without its seed**, which made it unreproducible (it is
phase D's seed-7 / seed-11 pair); and the C0-shape numbers **68** and
**≈ 260** were read at four sites as the *value* of σ^el_T/σ^el_U at x = 0.10
when they are its *suppression* against POLRAD's deuteron +0.064. It caught
one stale claim of a different kind as well: `PHYSICS_CHANNELS.md` still said
**no C++ path** opens `deuteron/fdeut.av18` or `li6.adr.fit`, and three paths
and one do, since C5.4 (§10 and §8).

## What the run did NOT establish

1. **No detection efficiency exists below Q² = 0.1 in any source this tree has
   seen**, and **no decay-lepton reconstruction efficiency exists in it at
   all** — only the geometry is bounded (A_geom = 0.994 at 10 × 99.5).
   **Still true after 2026-09-15**, when the sibling `../PolarizedLithiumSim`
   was surveyed for one — every module of `fastsim/polli_fastsim/`, all of
   `tools/fullsim/` and `tools/analysis/`, `plans/03`, `plans/09` and all 54
   `refs/` entries — and supplied none either: one per-track stand-in
   (`evgen/polligen/hfs.py:251`), one *constructed* scattered-electron ID
   profile (`evgen/polligen/reco.py:482`), no J/ψ reconstruction and no
   EICrecon output at all, so the factor stays **unbounded, not bounded**
   (`run_2026-09-06/phase_C_numbers.md` §C1). The
   far-forward working point is unchosen, and `Optics::lumi_fraction` is the
   single correction that on its own **restores a NO** (0.73 … 0.92 σ,
   106–167 fb⁻¹/u).
2. **Whether O5's band TOP crosses 3 σ.** Five defensible forms give 2.84,
   2.86, 3.03, 3.04, 3.29 σ; the crossing is withdrawn and the self-check now
   fails if it is re-asserted. The low edge and "3 σ inside {1, 10, 100}" hold.
3. **b₁ beyond A = 2 is not validated.** The gate is A = 2; the α–d step is
   untested, so the ±100 % band stays, and the *shipped* unpolarised default is
   outside the window the lift was granted on.
4. **b₁(⁷Li) is not implemented**, deliberately: its sign flips with the
   unpolarised backend. ~~and the Q(⁷Li) it would be gated against is quoted
   from memory and absent from this tree.~~ **The second half is CLOSED
   2026-09-06** (run 2026-09-06 task B3): Q(⁷Li) is sourced and in the tree as
   `LI7_QUADRUPOLE_FM2` = **−4.06 fm²** (`rc.hpp`, from the same TUNL A = 5,
   6, 7 evaluation as `LI6_QUADRUPOLE_FM2`), and `li7_alpha_t_quadrupole`
   makes the A = 7 gate a real, pinned test — **reported ratio 0.858389,
   14 % low** (0.871265 against the −4.00(3) previously quoted). The first
   half stands: that gate validates the α–t **wave function**, so **b₁(⁷Li)
   is still not implemented** and still blocked on the unpolarised-backend
   decision.
5. ~~**The exact Mo–Tsai elastic tail was never obtained.** What ships is a band
   of two stated models; the **tensor** fraction of the s-/p-peaks is *unknown,
   not zero*;~~ **REVISED 2026-09-06** (restored here 2026-09-15 as a
   strikethrough: this item had been rewritten in place, so the 2026-09-05
   reading was not recoverable from the tree — the convention, stated on run
   2026-09-06's own CW row, is that *dated records keep their originals and gain
   a correction box*, which item 4 above follows and this one did not. The
   pointer for what changed is
   `../run_2026-09-06/phase_B_numbers.md` §B1/§B2.) **What stands now:**

   **The exact Mo–Tsai elastic tail was never obtained** — and still has not
   been: `[MT69]` is not in this tree and no number from it is quoted anywhere.
   What ships is **three** tail models since 2026-09-06, when POLRAD's own
   exact one (`--rc-tail-model polrad-full`, Eq. (18) + Appendix B +
   Eq. (A.4)) was implemented; it is checked POLRAD-internally and against the
   leading-log fallback, **not** against an external exact tail, and it is
   **outside** the two t-peak models on 43.5 % of the sampler's cells, so
   those two are a **price range and not a confidence interval**
   (`../run_2026-09-06/phase_B_numbers.md` §B2). The **tensor** fraction of
   the s-/p-peaks is *computed* on that model and, on the two t-peak ones,
   since 2026-09-06, *bounded, not computed* — and the bound is **partial and, at the
   standard points, empty**: `--rc-sp-tensor-scale` (default 0) prices the
   **elastic** s/p column, which recovers **0.0 %** of the tensor-fraction
   collapse ×0.66139 / ×0.0031829 / ×6.6076e−05 at x = 0.01 / 0.10 / 0.30,
   Q² = 5 (bit-identical there: ⁶Li's coherent form factor is 45+ decades down
   at the s/p vertex Q′² ≈ Q²) and **at most 21.72 %** of it anywhere on the
   grid; the **quasi-elastic** s/p column, which is the whole of that collapse
   at x ≥ 0.10, is bounded by **nothing**
   (`../run_2026-09-06/phase_B_numbers.md` §B1) **and is not computed by the
   exact model either** — a nucleon has no tensor structure function, so the
   quasi-elastic tail is tensor-blind at all three peaks on all three models
   (§B2); the polarised quasi-elastic
   tail is *priced, not computed*; no RC
   calculation exists for a tagged tensor asymmetry, and everything the band is
   calibrated on is **deuteron**.
6. **O3 in the form asked** — the α–d split from GFMC configurations — remains
   unanswerable here: this tree has no GFMC configurations. What exists is a
   bound on the *size* of the effect.
6b. **A shipped tensor number was wrong and is now fixed — found after this
   run closed.** On 2026-09-06 the tagged sector's S–D interference sign was
   shown inverted: `TaggedModel::build_amp2` summed ψ_L where a partial-wave
   amplitude needs φ_L = i^L ψ_L, so the relative phase `(-1)^floor(L/2)` was
   missing (`../../benchmarking/07_cw_sign_investigation.md`, three
   independent derivations). **⁶Li A_zz^tag(k = 0.20 GeV), acceptance-weighted
   at the YR high-acceptance optics, moved +0.8450 → −1.2069 (Hulthén β = 0.30)
   and +0.4518 → −0.5191 (VMC AV18)** — the corrected VMC column is, digit for
   digit, the `vmc-flipD` *control* column this run published. Cosyn–Weiss II
   Eq. (6.12) now holds as an identity, `max|A_zz − CW| = 8.881784e−16` over
   26 880 cells against **2.740499e+00** before. **Nothing in this summary's
   own numbers moves**: ⁷Li does not move at all (one L = 1 wave), the
   spin-blind rate does not move (≤ 4.0e−16 relative), and every dilution,
   effective polarization and b₁ number here is angle-integrated or a closed
   form in P_D. Record: `../run_2026-09-06/phase_CW_numbers.md`.

7. **Nothing was pushed, published or sent.** `ci.yml` has never been parsed by
   GitHub (no cache round-trip, no runner-image apt set, no real 4-core
   timings); no wheel was uploaded, and its two data trees are still outside it
   with the build machine's absolute paths baked in; the Mäntysaari-group letter
   is drafted and unsent.
8. **Three costs were deliberately left unpriced** rather than estimated, and
   the registry says so in the word: the third R option for
   `Li6ConvolutionOptions` — one R hook threaded through both the numerator
   and the denominator, which no run was ever made with (row 3); the RC band's
   anchor alternatives 0.266 / 0.113, named with no band re-evaluation
   (row 17); and honouring `--pzz` on `helicity-flip`, which would move a
   shipped fill with no measured effect on any observable (row 20 / §15.5 D13).
   **ALL THREE ARE PRICED SINCE 2026-09-06** (run 2026-09-06 phase A,
   `cdd8591`): each was *implemented opt-in at a no-op default and run*, not
   estimated — `--r-source {unset,sigma-lt,r1998}` (row 3; §A1), the three-row
   δ_low band table (row 17; §A2, `python/tests/test_rc_low_x_anchor.py`) and
   `--pzz-mode {ladder,typed}` (row 20 / D13; §A3), all in
   `../run_2026-09-06/phase_A_numbers.md`. The word "unpriced" is gone from
   all three rows, from §B7(d), §B13(d) and §B18(d), and from
   `AUTHOR_DECISIONS.md`'s Part B item 8. **None of the three was decided**:
   pricing a cost is not taking it, and all three rows are still open —
   row 3 has in fact become a *narrower* question, *(iii) with which R?*
9. **Twenty-five decisions sat in the registry and this phase took none of
   them.** Six were marked *applied* by the phase that applied them — four
   *applied — confirm or revert* (rows 1, 11, 18, 19) and two *applied —
   confirm* (rows 4, 12); the rest were open. All are listed in `STATUS.md`'s
   table and stated in full in `AUTHOR_DECISIONS.md`.
   **The registry is TWENTY-SEVEN rows since 2026-09-15**, and the follow-on
   run took none of them either. Row **26** (§B27) is the tagged S–D
   interference phase and the `tagged.json` re-pin, *applied — confirm or
   revert*; row **27** (§B28) is the sibling repository's published
   `A_zz^tag` numbers, open. Row **6** also became *applied in part — confirm*
   on 2026-09-06, when the corrected `A_zz(Born)` column was republished
   beside the original. Counted at HEAD from the index table: **eight rows are
   marked applied** — five *confirm or revert* (1, 11, 18, 19, 26), two
   *confirm* (4, 12) and one *applied in part — confirm* (6) — one is
   *decided on the record* (7), three are *carried forward, confirm only*
   (22, 23, 24) and **fifteen are open** (2, 3, 5, 8, 9, 10, 13, 14, 15, 16,
   17, 20, 21, 25, 27).
