# Phase C — the measured numbers of `design_C_tensor_rc.md` §8

Everything below is a **measurement of the shipped model**, not a re-statement
of the design's estimates. Where the two disagree, the measurement replaces the
estimate and the row says so.

> **Revision of 2026-09-03 (physics + integration review).** Every number in
> §8.1, §8.1c and §8.3 below was **regenerated**. Two things moved them:
> (i) the **per-nucleon reduction of the elastic tail was wrong by a factor
> `A = 6`** — Eq. (38) is the *whole-nucleus* `d²σ/(dx_A dy)` and needs
> `(1/A)·(dx_A/dx) = 1/A²`, not the single `m_p/M_A` v0 applied (see
> `polrad_transcription_check.md` §4, corrected); and (ii) the **quasi-elastic
> Pauli suppression `S(q)` is now on by default** (`RcOptions::qe_kf_gev =
> 0.169 GeV`, POLRAD Eq. (44) / `ffquas`), which was design Q9 and is now
> closed. The nuclear map also moved from POLRAD's fixed-target
> `x_A = x m_p/M_A` to this library's exact `x_A = x/A` (0.3–7 % on `σ^el_U`).
> **Everything the old file said about the tail's SIZE is superseded.**

**Configuration for every table.** ⁶Li, `--config 1` (10 GeV e × 99.5 GeV/u,
`s_per_nucleon = 3980.0000 GeV²`), `--channel inclusive`,
`--plan tensor-thirds --pzz 0.6`, `--rc tensor-band` at the defaults,
`θ_S = 0`, `RcTailModel::TPeak`, `n_eta = 128`, the production `100 × 72`
sampler grid on `generator_scenario` (`Q² ≥ 0.7`, `y ∈ [0.004, 0.985]`,
`W² ≥ 8`), 3051 accepted cells. The `b₁` input is `toy_b1` (the digitised
default), so every `A_zz` below is *this generator's model* `A_zz` and not a
measurement — the RC quantities are what this file is about.

Tail node grid: **101 × 77** in `(ln x, ln y)` — the sampler's own `x` cell
**edges** and a geometric `y` grid over the scenario window, with the design's
refinement rows at `y ∈ {0.9, 0.95, 0.97, 0.98, 0.985}` merged in.

---

## 8.1 The radiative tails

### 8.1a `w_tail − 1` (elastic + unpolarised quasi-elastic, as a fraction of the Born)

`w_tail − 1 = [σ^el_U + (Q_N/6)σ^el_T + κ_qe σ^q_U] / [x·s·dσ_unpol·(1 + w_avg)]`
with `Q_N = P_zz^eff` (+1 at `m = ±1`, −2 at `m = 0`) and `κ_qe = 1`.

| x | Q² = 2 GeV² (y) | Q² = 5 GeV² (y) | Q² = 10 GeV² (y) |
|---|---|---|---|
| 0.01, `m = ±1` | **8.9654e−05** (0.0503) | **5.18398e−04** (0.1256) | **2.14577e−03** (0.2513) |
| 0.01, `m = 0`  | **8.95284e−05** | **5.1787e−04** | **2.14413e−03** |
| 0.10, `m = ±1` | **2.18453e−07** (0.0050) | **1.32298e−06** (0.0126) | **5.19642e−06** (0.0251) |
| 0.10, `m = 0`  | **2.1674e−07** | **1.31299e−06** | **5.15836e−06** |
| 0.30, `m = ±1` | *does not exist* (y = 0.0017 < 0.004) | **6.8624e−08** (0.0042) | **2.72302e−07** (0.0084) |
| 0.30, `m = 0`  | — | **6.8605e−08** | **2.72227e−07** |

> **The headline, and it is not what a fixed-target intuition expects.** At EIC
> collider kinematics these `(x, Q²)` cells sit at **small y** (`y = Q²/(x s)`
> with `s = 3980`), and the tail carries `Y₊ = [1 + (1−y)²]/(1−y) ≈ 2` there.
> The radiative-tail dilution of `A_zz` is therefore **five to eight orders of
> magnitude smaller than at HERMES**, where the same `x` is reached at
> `y ≈ 0.5`. §8.1b is the same model at HERMES-like `y` and gives 1.8 % at
> `y = 0.5` and 33 % at `y = 0.9`. **`rc_tail` is a small effect at the EIC and
> a large one at a fixed target, and the difference is entirely `Y₊` and the
> Born's `1/Q⁴`.**
>
> **And it is a LOWER BOUND, because this is the t-peak alone.** `T8(c)` now
> measures the leading-log s- and p-peaks of the same observable
> (`tests/test_rc.cpp`, an independent Weizsäcker–Williams × elastic
> construction). For ⁶Li at EIC `Q² ≥ 20 GeV²` the t-peak is **> 99 %** of
> `t + s + p` for both tails — POLRAD §2.1.3 B's claim, confirmed *in the
> regime this generator runs in*. It is **not** true elsewhere: at the HERMES
> deuteron point (`x = 0.012`, `y = 0.85`, `Q² = 0.53 GeV²`) the t-peak is only
> **23 %** of the leading-log total (elastic 65 %, quasi-elastic 17 %), i.e.
> low by a factor **4.4** — which is what reconciles the ~1 % the t-peak alone
> gives there with HERMES's "almost 50 % of the statistics in the lowest-x
> bin". At the low-`Q²` corner of the generator window (`x = 0.01`, `y = 0.1`,
> `Q² = 4 GeV²`) the quasi-elastic s-peak is already **3.3×** the t-peak.
> **Never quote `rc_tail` as "the" radiative tail.**

### 8.1b The `y`-scan at fixed `Q² = 5 GeV²` (the scan T5 is written against)

| y | 0.1 | 0.3 | 0.5 | 0.7 | 0.9 | 0.97 |
|---|---|---|---|---|---|---|
| x (= Q²/(y s)) | 0.012563 | 0.0041876 | 0.0025126 | 0.0017947 | 0.0013959 | 0.0012951 |
| `w_tail − 1` | **3.00722e−04** | **4.3188e−03** | **1.79688e−02** | **6.23201e−02** | **3.27277e−01** | **1.29791** |
| `Y₊(y)` | 2.011 | 2.129 | 2.500 | 3.633 | 10.100 | 33.363 |
| `(w_tail−1)/Y₊` | 1.4953e−04 | 2.02897e−03 | 7.18751e−03 | 1.71523e−02 | 3.24036e−02 | 3.89024e−02 |

**Design vs code, stated plainly.** §5's T5(i) asks that `(w_tail−1)/Y₊` "varies
by less than 3× across the scan". **It does not, and a correct implementation
cannot**: at fixed `Q²` the scan also sweeps `x` by a factor 10, and the
design's own §8.1 says "the x-dependence is the strong one". The measured
spread across the full scan is **292×**; across the **last three points** it is
**2.3×**, and the quantity is **monotone increasing in y** at every step.
`tests/test_rc.cpp` asserts the monotonicity everywhere and the "< 3×" only
where it is true, and says why in place.

### 8.1c The derived quantities that actually matter

| x | Q² | `(1/6)σ^el_T/σ^el_U` | `σ^q_U/σ^el_U` (S(q) on) | `σ^q_U/σ^el_U` (S = 1) | `A_zz(Born)` | `ΔA_zz` (tail) |
|---|---|---|---|---|---|---|
| 0.01 | 2 | **−5.0943e−04** | **0.28237** | 0.5981 | −1.7299e−03 | **+8.3816e−08** |
| 0.01 | 5 | −5.0943e−04 | 0.28237 | 0.5981 | −1.4735e−03 | **+3.5167e−07** |
| 0.01 | 10 | −5.0943e−04 | 0.28237 | 0.5981 | −1.3055e−03 | **+1.0939e−06** |
| 0.10 | 2 | **+1.5507e−04** | **2.7479** | 3.164 | −5.1730e−03 | +1.1451e−09 |
| 0.10 | 5 | +1.5595e−04 | 2.7475 | 3.164 | −4.9751e−03 | +6.6751e−09 |
| 0.10 | 10 | +1.5622e−04 | 2.7474 | 3.164 | −4.8235e−03 | +2.5436e−08 |
| 0.30 | 5 | **+1.3287e−02** | **1000.2** | 1000.2 | −1.5852e−04 | +1.2665e−11 |
| 0.30 | 10 | +1.3464e−02 | 975.3 | 975.3 | −1.5790e−04 | +5.0374e−11 |

`ΔA_zz = (A_zz + 2 r_T)/(1 + r_U) − A_zz` with `r_U` the unpolarised tail ratio
and `r_T` its `Q_N`-coefficient — i.e. the shift a fit that ignores the tail
would make, **dominated by the unpolarised dilution** exactly as the design
predicted.

Three things this table settles:

1. **The tensor fraction of the elastic tail changes SIGN with x** — −5.1e−04 at
   `x = 0.01`, +1.7e−04 at `x = 0.10`, +1.3e−02 at `x = 0.30`. This is the ⁶Li
   counterpart of the deuteron sign change the transcription check measured with
   POLRAD's own code (+0.106, +0.064, −0.117 at three points) and it is the
   reason **`σ^el_T` is a separate table and never a scale factor on `σ^el_U`**.
   Design §2.1 expected `O(10⁻²)` per unit `Q_N`; that is met at `x = 0.30` and
   is **20–60× too large at `x ≤ 0.1`**.
2. **`σ^q_U/σ^el_U` is 0.28 → 2.75 → 1000 over `x = 0.01 → 0.30`** (0.60 → 3.16
   → 1000 with the Pauli suppression off). Design §1.4.3 estimated "30–70 % of
   the ERT over x = 0.01–0.3"; **this table replaces that estimate**, and the
   estimate is right only at the very bottom of the range. **`rc_tail` is
   quasi-elastic-dominated at every `x` ≳ 0.03**: the QRT is 22 % of the tail at
   `x = 0.01`, 73 % at `x = 0.1` and **99.9 %** at `x = 0.30`. So
   `--rc-qe-suppression` and `RcOptions::qe_kf_gev` are the dominant tail knobs
   almost everywhere, and `qe_suppression = 0` is never a small variation. (The
   pre-2026-09-03 numbers here read 0.10 / 0.53 / 178 — six times too small,
   because the elastic tail carried one per-nucleon factor too few.)
3. **The tail is not the leading RC systematic at the EIC.** `ΔA_zz` from the
   whole tail is `3.5e−07` at `x = 0.01, Q² = 5`, against a **band half-width of
   `4.4e−04`** (§8.3) — a factor 1300. The band, i.e. the *unapplied
   lepton-vertex correction*, is what `--rc tensor-band` is for; the tail is a
   bookkeeping item at these energies and the headline at a fixed target. Note
   the tail's own factor-4 s-/p-peak deficit at fixed-target kinematics (§8.1a's
   second box) does **not** change that ordering here.

---

## 8.2 The band — `w_hi − 1` at `P_zz = +1`, `θ_S = 0`, `Q² = 5 GeV²`

`τ = (P_zz/2)A_zz / [1 + (P_zz/2)A_zz]` — **with** the denominator the first
draft omitted.

| x | δ(x) | `A_zz` | `τ` | `w_hi − 1` | `w_lo − 1` |
|---|---|---|---|---|---|
| 0.010 | 0.30000 | −1.4735e−03 | −7.3728e−04 | **−2.2119e−04** | +2.2119e−04 |
| 0.063 | 0.11081 | −4.2848e−03 | −2.1470e−03 | **−2.3790e−04** | +2.3790e−04 |
| 0.100 | 0.06331 | −4.9751e−03 | −2.4937e−03 | **−1.5789e−04** | +1.5789e−04 |
| 0.160 | 0.01500 | −4.1666e−03 | −2.0876e−03 | **−3.1315e−05** | +3.1315e−05 |
| 0.300 | 0.01500 | −1.5852e−04 | −7.9268e−05 | **−1.1890e−06** | +1.1890e−06 |

`w_lo + w_hi == 2.0` bit-for-bit on every event (T4a). The band peaks near
`x = 0.063` — not at the lowest x — because `δ(x)` is still rising there while
`|A_zz|` has not yet fallen: **the largest RC systematic on `A_zz` sits in the
middle of the low-x range, not at its edge.**

---

## 8.3 Run-level summary

| quantity | value |
|---|---|
| `RcModel` construction (101 × 77 nodes, `n_eta = 128`, 4 quadratures/node) | **0.106 s** — the `Pipeline` setup difference `--rc tensor-band` minus `--rc off`, median of 7 (0.1416 s vs 0.0353 s). **This is the one construction-time number**; `OPEN_ITEMS_SOLUTIONS.md` §9 and `USAGE.md` §7b quote the same measurement, and the earlier 0.134 s / 0.16 s pair is withdrawn. |
| `RcModel::fill` cost per event (slot 0 only, 3051-cell loop × 200) | **0.206 µs/event = 4.86e6 ev/s** |
| throughput with `--rc tensor-band` vs `--rc off`, END TO END | **574 852 → 495 017 ev/s = −13.9 %** (`Pipeline.generate`, 400 k inclusive events, single thread, best of 3, same build and machine as the construction row). The absolute rates are machine- and load-dependent; the **−14 %** is the number to quote. |
| `clipped_cell_fraction()` at `tail_max = 10`, global | **0 %** (was 1.71 % before the per-nucleon fix, which lowered the whole elastic tail by A = 6) |
| `clipped_fraction_by_y()` — `y < 0.5` / `0.5–0.9` / `> 0.9` | **0 % / 0 % / 0 %** |
| clipped EVENTS, default 2000-event inclusive CLI run — tail / band | **0 / 0** (`meta["rc_clipped_tail_event_fraction"]`, `..._band_...`; before the fix one event in 2000 hit `tail_max`). This is a **different quantity** from the node fractions above: nodes are not event-weighted, and the node statistic is computed at `q_n = 0` while the per-event clip carries the `(q_n/6)` tensor term. |
| clipped EVENTS on the BAND, 20 k tagged-alpha (`--pzz 0.6`) | **0.62 %** on ⁶Li, **2.6 %** on ⁷Li (API path). `tau_tag = 1 − n̄/n_M` reaches **30.7** at the nodes of the M = 0 spectator density; `RcOptions::band_tau_max = 1` is what keeps every published edge inside `[1 − δ, 1 + δ] = [0.7, 1.3]` instead of the **−1.79** (⁶Li) / **−8.35** (⁷Li) the unclamped v0 emitted. |
| fitted `(a, α, q₀)` and χ²/ndf against the Li/Suelzle data | **NOT DONE.** `a = 1.9069 fm`, `α = 0.13822`, `q₀ = 3.0999 fm⁻¹` are the design's **unfitted starting values**, and they are shipped as such: the Suelzle–Yearian–Crannell (PR 162 (1967) 992) and Li–Sick–Whitney–Yearian (NPA 162 (1971) 583) elastic data are **not in this repository in any machine-readable form**, so there is nothing to fit against. They do reproduce their two design targets exactly (`⟨r²⟩_point = 6.0788 fm²`, first C0 zero at 3.0999 fm⁻¹, both asserted in T11). **Q1 stays open.** |
| fitted `(q_z, b)` of `F_mag` against WS98's three `F_T` landmarks | **NOT DONE.** `q_z = 1.30 fm⁻¹`, `b = 1.85 fm`. WS98 (nucl-th/9807037) is **figures only — there is no table anywhere in the paper** (design §2.1), and a two-parameter shape cannot hit three landmarks in any case. `--rc-tail-tensor-scale 0.5 / 1 / 2` bands the NORMALISATION, not the shape. **Q10 stays open.** |
| δ(A_zz) at `x = 0.01, Q² = 5` from the band, `δ_low = 0.19` vs `0.30` | **2.7996e−04 vs 4.4204e−04** (half-widths on `A_zz`) — **the dominant RC systematic, by a factor 300 over everything below** |
| δ(A_zz) at `x = 0.01, Q² = 5` from `fq_scale = 0` vs `2` | `ΔA_zz` **+7.54154e−07 → −4.9925e−08**, spread **8.041e−07** |
| δ(A_zz) at `x = 0.01, Q² = 5` from `tail_tensor_scale = 0.5` vs `2` | `ΔA_zz` **+3.59193e−07 → +3.23118e−07**, spread **3.61e−08** |
| δ(A_zz) at `x = 0.01, Q² = 5` from `qe_suppression = 0` vs `1` | `ΔA_zz` **+1.83644e−07 → +3.51674e−07**, spread **1.680e−07** |
| δ(A_zz) at `x = 0.01, Q² = 5` from `qe_kf_gev = 0` (S = 1) vs `0.169` | `ΔA_zz` **+5.39491e−07 → +3.51674e−07**, spread **1.878e−07** — the Pauli suppression is a *larger* knob than `qe_suppression = 0 → 1` at this x, which is why it is the default rather than a band edge |

**The `fq_scale` band CHANGES THE SIGN of `ΔA_zz`** (+3.9e−06 at 0, −9.6e−07 at
2) and is **not** symmetric about the nominal — `σ^el_T` is *quadratic* in
`F_q`. That is exactly why the design forbids rescaling one run and why T12 fits
a quadratic through three runs instead of asserting linearity.

**Clipping.** **Nothing clips any more.** Before the 2026-09-03 per-nucleon fix
21.95 % of the `y > 0.9` nodes exceeded `tail_max = 10` and one event in 2000
hit the ceiling in the default CLI run; dividing the elastic tail by the missing
`A = 6` removed all of it. The `y > 0.9` corner is still where `Y₊ ~ 1/(1−y)`
would bite first, so the four node fractions and the two **event** fractions
(`meta["rc_clipped_tail_event_fraction"]`, `..._band_event_fraction"]`, and
`Event.rc_clipped` per event) are all still printed by the run. What *does*
clip now is the **band on the tagged channels**, where `tau_tag` diverges at the
`n_M` nodes — 0.62 % of ⁶Li tagged-alpha events, 2.6 % of ⁷Li.

**Form factor in use** (printed by the run banner):

```
HoSpin1FF (design_C_tensor_rc.md sec. 2.1, PHENOMENOLOGICAL FIT, not a shell
model): F_point(a = 1.9069 fm, alpha = 0.13822), F_mag(q_z = 1.3 fm^-1,
b = 1.85 fm); MEASURED moments mu = 0.822047 mu_N, Q = -0.0818 fm^2 (NOT VMC:
WS98's Q(6Li) = -0.23(9) fm^2 is 3x the measured one); Z = 3,
M_A = 5.60152 GeV (beams.cpp AME mass); folded with G_E^p + G_E^n (N = Z
isoscalar); fq_scale = 1, tail_tensor_scale = 1.  UNFITTED STARTING VALUES ...
```

---

## Honest flags to repeat wherever any number above is quoted

* `RC_DELTA_LOW_X = 0.30` is the **size of a correction this generator does not
  apply**, taken as a 1σ band. Its source (Gakh–Shekhovtsova, hep-ph/0403262)
  has **zero INSPIRE citations**. Band it: run 0.19 and 0.30, never quote one.
* **Everything cited is DEUTERON.** Whether the deuteron's *fractional* RC
  transfers to ⁶Li is untested and is the largest unquantified assumption here.
* The ⁶Li form-factor parameters are **unfitted starting values**, not a fit
  (two rows of §8.3). `F_q`'s shape between its two pinned ends is model, and
  `σ^el_T` is quadratic in it.
* **No RC calculation exists for a tagged tensor asymmetry**, and none for any
  φ-dependent tensor observable at any axis.
* The **polarised** quasi-elastic tail is not priced at all (Zhou et al.,
  PRL 82 (1999) 687), and at `x = 0.30` the tail is **99.9 %** quasi-elastic.
* **`rc_tail` is the t-PEAK ONLY and is therefore a LOWER BOUND.** Measured
  against the leading-log s-/p-peaks (T8(c)): fine for ⁶Li at `Q² ≥ 20 GeV²`
  (> 99 %), **low by a factor 4.4 at the HERMES deuteron point** and by 3.3× on
  the quasi-elastic piece at the `Q² ≈ 4 GeV²` corner of the generator window.
* The **tagged band is clamped** (`RcOptions::band_tau_max = 1`). `tau_tag`
  itself is unbounded at the nodes of `n_M(k, c)`; the clamp is a *choice*, the
  clipped fraction is reported, and `n_M → 0` is exactly where the
  fractional-rescale ansatz stops meaning anything.
