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
> **And the DEFAULT is a LOWER BOUND, because it is the t-peak alone — with
> the other edge shipped since the 2026-09-03 run.** `T8(c)` measures the
> leading-log s- and p-peaks of the same observable (an independent
> Weizsäcker–Williams × elastic construction), and since B2 it measures them
> through the **shipped** `ll_peaks_spin1` / `ll_peaks_qe` rather than a
> test-local copy; `RcTailModel::TPeakPlusLL` (`--rc-tail-model t-peak+ll`)
> puts them in the tail. **The two are a band.** The upper edge is a **STATED
> MODEL of mixed approximation orders** — an η_A quadrature plus a single-z
> collinear leading log, good to the worse of the two (~5–10 %), with an
> uncancelled soft `1/(1−z)` as `y → 0` — not a controlled O(α) expansion; and
> **neither edge has a tensor s/p peak**, so `TPeakPlusLL` **lowers the tensor
> fraction** of the tail (that fraction is unknown, not zero).
>
> POLRAD §2.1.3 B's claim is an **event-weighted** statement about the bulk of
> the regime this generator runs in and is **false cell by cell** (corrected
> 2026-09-04; `phase_B_numbers.md` §B2.3(i)). For ⁶Li at EIC `Q² ≥ 20 GeV²`
> and `y ≤ 0.9` the *event-weighted* mean dilution moves only **+0.61 %**
> (8.48914e−03 → 8.54072e−03), but **331 of that window's 1356 accepted cells
> — 24.4 %, 28.2 % of its cross section — differ by more than 1 %**, worst
> **×6444** at `x = 0.7943`, `y = 0.0088`, and **318 of them sit at `y < 0.1`**.
> Per cell the two edges agree to **0.55 %** only for `0.15 ≤ y ≤ 0.7`
> (0.77 % to `y ≤ 0.8`, 2.0 % to `y ≤ 0.9`); below `y = 0.15` there is no
> agreement statement. Two corners, two mechanisms: at `y = 0.985` they differ
> by **59 %** because `z_s = (1−y)/(1−x_A y) → 0` puts the elastic vertex at
> `Q′² → 0` and the s-peak beats `Y₊`; at `y → 0` it is the **uncancelled soft
> radiator**, `1 − z_s = y(1−x_A)/(1−x_A y)` so `D(z_s) ∝ 1/y` (11.5 at
> `x = 0.72` on the `Q² = 23.88` line), against a t-peak whose own vertex is
> stuck above `t_min ∝ x_A²` and which the `Q² ≥ 20` cut pushes to high `x`.
> That corner is a **breakdown of the upper edge**, not evidence that the
> t-peak is low by ×6444. It is **not** true elsewhere either: at the HERMES
> deuteron point (`x = 0.012`, `y = 0.85`, `Q² = 0.53 GeV²`) the t-peak is only
> **23 %** of the leading-log total (elastic 65 %, quasi-elastic 17 %), i.e.
> low by a factor **4.36** — which is what reconciles the ~1 % the t-peak alone
> gives there with HERMES's "almost 50 % of the statistics in the lowest-x
> bin". At the low-`Q²` corner of the generator window (`x = 0.01`, `y = 0.1`,
> `Q² = 4 GeV²`) the quasi-elastic **s+p** is **3.35×** the t-peak
> unsuppressed and **7.09×** at the shipped Pauli `k_F`; the s-peak *alone* is
> **2.33×** / **4.94×**, and the earlier "s-peak … 3.3×" wording here conflated
> the two. **Never quote either edge as "the" radiative tail.**

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

> **BANDED 2026-09-04 (task B1).** Every `σ^el_T` row below was made with the
> **unfitted** `C0Shape::Ho` monopole. The second edge of that shape band
> (`--rc-c0-shape vmc-ft`, the j₀ transform of `data/vmc/density/li6.density` at
> the *same* measured ⟨r²⟩_point) is added as a second row at every point.
> **The `Ho` numbers are unchanged and reproduce exactly**; what changed is what
> may be concluded from them — see the rewritten item 1 and
> `../run_2026-09-03/phase_B_numbers.md` §B1.

| x | Q² | edge | `(1/6)σ^el_T/σ^el_U` | `σ^q_U/σ^el_U` (S(q) on) | `σ^q_U/σ^el_U` (S = 1) | `A_zz(Born)` | `ΔA_zz` (tail) |
|---|---|---|---|---|---|---|---|
| 0.01 | 2 | ho | **−5.0943e−04** | **0.28237** | 0.5981 | −1.7299e−03 | **+8.3816e−08** |
| 0.01 | 2 | vmc-ft | **−5.4705e−04** | 0.27942 | 0.5918 | −1.7299e−03 | +7.9026e−08 |
| 0.01 | 5 | ho | −5.0943e−04 | 0.28237 | 0.5981 | −1.4735e−03 | **+3.5167e−07** |
| 0.01 | 5 | vmc-ft | −5.4705e−04 | 0.27942 | 0.5918 | −1.4735e−03 | +3.2291e−07 |
| 0.01 | 10 | ho | −5.0943e−04 | 0.28237 | 0.5981 | −1.3055e−03 | **+1.0939e−06** |
| 0.01 | 10 | vmc-ft | −5.4705e−04 | 0.27942 | 0.5918 | −1.3055e−03 | +9.7195e−07 |
| 0.03 | 5 | ho | −7.3583e−04 | 0.54628 | 0.8679 | — | — |
| 0.03 | 5 | vmc-ft | −8.0426e−04 | 0.53452 | 0.8492 | — | — |
| 0.10 | 2 | ho | **+1.5507e−04** | **2.7479** | 3.164 | −5.1730e−03 | +1.1451e−09 |
| 0.10 | 2 | vmc-ft | **−4.1523e−05** | 2.4884 | 2.865 | −5.1730e−03 | +1.1531e−09 |
| 0.10 | 5 | ho | +1.5595e−04 | 2.7475 | 3.164 | −4.9751e−03 | +6.6751e−09 |
| 0.10 | 5 | vmc-ft | **−4.0702e−05** | 2.4881 | 2.865 | −4.9751e−03 | +6.7163e−09 |
| 0.10 | 10 | ho | +1.5622e−04 | 2.7474 | 3.164 | −4.8235e−03 | +2.5436e−08 |
| 0.10 | 10 | vmc-ft | **−4.0458e−05** | 2.4880 | 2.865 | −4.8235e−03 | +2.5575e−08 |
| 0.30 | 5 | ho | **+1.3287e−02** | **1000.2** | 1000.2 | −1.5852e−04 | +1.2665e−11 |
| 0.30 | 5 | vmc-ft | +1.5556e−02 | 329.85 | 329.85 | −1.5852e−04 | +1.7334e−11 |
| 0.30 | 10 | ho | +1.3464e−02 | 975.3 | 975.3 | −1.5790e−04 | +5.0374e−11 |
| 0.30 | 10 | vmc-ft | +1.5622e−02 | 323.66 | 323.66 | −1.5790e−04 | +6.9223e−11 |


> **⚠ CORRECTION 2026-09-04 (task B3): the `A_zz(Born)` column of this table and
> of §8.2 is ×3.253983 too large, and the `ΔA_zz` column moves with it — at
> x = 0.01 it changes SIGN.** That column reproduces exactly as
> `azz(toy_b1(x, q2, f1, B1Mode::Digitized), …)` — the **deuteron** Miller b₁
> table — while the shipped `default_inclusive_kernel` has used
> `Li6B1(MillerB1)` for every spin-1 ion since before the RC commit. Measured
> at Q² = 5: this generator's own ⁶Li `A_zz(Born)` is
> **−4.528242e−04 / −1.528923e−03 / −4.871661e−05** at x = 0.01 / 0.10 / 0.30,
> against the −1.4735e−03 / −4.9751e−03 / −1.5852e−04 printed here. So the
> §8.2 band half-width at x = 0.01 is **1.3585e−04**, not 4.4204e−04, and
> `ΔA_zz` at the defaults and x = 0.01 is **−1.76980e−07**, not +3.51674e−07
> — the sign flips because ΔA_zz = [2 r_T − A_zz r_U]/(1 + r_U) and the two
> terms are comparable there. The x = 0.10 and x = 0.30 rows keep their sign
> and shrink (**+2.127452e−09** and **+5.130156e−12**, measured, against
> +6.6751e−09 and +1.2665e−11 printed). **The `(1/6)σ^el_T/σ^el_U`,
> `σ^q_U/σ^el_U` and `w_tail` columns
> are NOT affected** and reproduce here to every printed digit. This table was
> deliberately **not** republished by task B3 — the finding is recorded with its
> cause so that regenerating it is a deliberate step.
> `../run_2026-09-03/phase_B_numbers.md` §B3.2.
>
> **REPUBLISHED BESIDE THE ORIGINAL ON 2026-09-06 (registry row 6 / §B11).**
> The table above **stays**, labelled as what it was computed with; the
> recomputation is the separate table below. Nothing was deleted and no
> decision was taken: option (i) "leave with the note" and option (ii)
> "republish" both remain open, and (ii) is now costed at zero because the
> corrected column exists beside the original.

The `x = 0.03` row is new (it is where the `Ho` edge's sign change happens) and
carries no `A_zz(Born)`, so it carries no `ΔA_zz` either. `σ^q_U` itself is
**bit-identical** between the two edges at every point — this knob is the
elastic C0 shape and nothing else — so the `σ^q_U/σ^el_U` column moves only
through its denominator, which rises by ×1.011 / ×1.022 / ×1.104 / ×3.03 at
x = 0.01 / 0.03 / 0.10 / 0.30.

`ΔA_zz = (A_zz + 2 r_T)/(1 + r_U) − A_zz` with `r_U` the unpolarised tail ratio
and `r_T` its `Q_N`-coefficient — i.e. the shift a fit that ignores the tail
would make, **dominated by the unpolarised dilution** exactly as the design
predicted.

Three things this table settles:

1. **The tensor fraction of the elastic tail changes SIGN with x — but WHERE it
   changes is not determined, and at `x = 0.10` the sign itself is a band edge.**
   Both edges of the C0 shape band are negative at `x = 0.01` and `0.03` and
   positive at `x = 0.30`, so the sign *change* is not an artefact of the
   unfitted shape; but the crossing sits between `x = 0.03` and `0.10` on `Ho`
   and between `0.10` and `0.30` on `vmc-ft`, and at `x = 0.10` the two edges
   read **+1.5595e−04 and −4.0702e−05** at every `Q²` in the table. **The
   +1.5595e−04 that this document published on 2026-09-02 may not be quoted as a
   result.** Adding the `fq_scale` band widens x = 0.10 to −9.0e−05 … +4.0e−04,
   still spanning zero on both edges. This is the ⁶Li counterpart of the
   deuteron sign change the transcription check measured with POLRAD's own code
   (+0.106, +0.064, −0.117 at three points) and it is the reason **`σ^el_T` is a
   separate table and never a scale factor on `σ^el_U`**. Design §2.1 expected
   `O(10⁻²)` per unit `Q_N`; that is met at `x = 0.30` and is **20–60× too large
   at `x ≤ 0.1`** — on **both** edges, so that part of the conclusion stands.
   Full band, its provenance and its two policy overrides:
   `../run_2026-09-03/phase_B_numbers.md` §B1.
2. **`σ^q_U/σ^el_U` is 0.28 → 2.75 → 1000 over `x = 0.01 → 0.30`** (0.60 → 3.16
   → 1000 with the Pauli suppression off; 0.28 → 2.49 → 330 on the `vmc-ft`
   edge, which does not change the reading). Design §1.4.3 estimated "30–70 % of
   the ERT over x = 0.01–0.3"; **this table replaces that estimate**, and the
   estimate is right only at the very bottom of the range. **`rc_tail` is
   quasi-elastic-dominated at every `x` ≳ 0.03**: the QRT is 22 % of the tail at
   `x = 0.01`, 73 % at `x = 0.1` and **99.9 %** at `x = 0.30`. So
   `--rc-qe-suppression` and `RcOptions::qe_kf_gev` are the dominant tail knobs
   almost everywhere, and `qe_suppression = 0` is never a small variation. (The
   pre-2026-09-03 numbers here read 0.10 / 0.53 / 178 — six times too small,
   because the elastic tail carried one per-nucleon factor too few.)
3. **The tail is not the leading RC systematic at the EIC.** `ΔA_zz` from the
   whole tail is `3.5e−07` at `x = 0.01, Q² = 5` (`3.2e−07` on the `vmc-ft`
   edge — the C0 band moves it by −8 % at `x = 0.01`, +0.6 % at `0.10` and
   +37 % at `0.30`, and changes nothing here), against a **band half-width of
   `4.4e−04`** (§8.3) — a factor 1300. **CORRECTED 2026-09-04 (task B6): the
   factor is 768.** Both entries used the ×3.253983 deuteron-b₁ `A_zz` of
   §8.2's warning; with the ⁶Li b₁ the tail is −1.769797e−07 and the band
   half-width 1.358472e−04 (×858 and ×1086 at `a_transfer_frac` 0.5 and 1).
   The ordering is unchanged. The band, i.e. the *unapplied
   lepton-vertex correction*, is what `--rc tensor-band` is for; the tail is a
   bookkeeping item at these energies and the headline at a fixed target. Note
   the tail's own factor-4 s-/p-peak deficit at fixed-target kinematics (§8.1a's
   second box) does **not** change that ordering here.

### 8.1c-corr `A_zz(Born)` and `ΔA_zz` recomputed with the SHIPPED kernel — 2026-09-06

**A dated correction, not a replacement.** §8.1c above is left exactly as it
was published on 2026-09-02/04. This table is the same eight `(x, Q²)` points
and the same two C0 edges with the **shipped** ⁶Li b₁ —
`default_inclusive_kernel`'s `Li6B1(MillerB1)`, i.e. `InclusiveKernel::tables`'
own `b1`/`b2` — in place of the `azz(toy_b1(x, q2, f1, B1Mode::Digitized), …)`
the original column reproduces. Only the two rightmost columns of §8.1c can
move: `(1/6)σ^el_T/σ^el_U`, `σ^q_U/σ^el_U` and `w_tail` carry no `A_zz` and are
unchanged (re-measured, identical to every printed digit).

| x | Q² | edge | `A_zz(Born)` **shipped ⁶Li b₁** | `ΔA_zz` **shipped ⁶Li b₁** | `A_zz(Born)` as published | `ΔA_zz` as published |
|---|---|---|---|---|---|---|
| 0.01 | 2 | ho | **−5.316313e−04** | **−2.355530e−08** | −1.7299e−03 | +8.3816e−08 |
| 0.01 | 2 | vmc-ft | **−5.316313e−04** | **−2.922677e−08** | −1.7299e−03 | +7.9026e−08 |
| 0.01 | 5 | ho | **−4.528242e−04** | **−1.769797e−07** | −1.4735e−03 | +3.5167e−07 |
| 0.01 | 5 | vmc-ft | **−4.528242e−04** | **−2.100985e−07** | −1.4735e−03 | +3.2291e−07 |
| 0.01 | 10 | ho | **−4.012087e−04** | **−8.419210e−07** | −1.3055e−03 | +1.0939e−06 |
| 0.01 | 10 | vmc-ft | **−4.012087e−04** | **−9.796923e−07** | −1.3055e−03 | +9.7195e−07 |
| 0.10 | 2 | ho | **−1.589729e−03** | +3.643984e−10 | −5.1730e−03 | +1.1451e−09 |
| 0.10 | 2 | vmc-ft | **−1.589729e−03** | +3.506747e−10 | −5.1730e−03 | +1.1531e−09 |
| 0.10 | 5 | ho | **−1.528923e−03** | +2.127452e−09 | −4.9751e−03 | +6.6751e−09 |
| 0.10 | 5 | vmc-ft | **−1.528923e−03** | +2.042099e−09 | −4.9751e−03 | +6.7163e−09 |
| 0.10 | 10 | ho | **−1.482347e−03** | +8.116142e−09 | −4.8235e−03 | +2.5436e−08 |
| 0.10 | 10 | vmc-ft | **−1.482347e−03** | +7.774146e−09 | −4.8235e−03 | +2.5575e−08 |
| 0.30 | 5 | ho | **−4.871661e−05** | +5.130156e−12 | −1.5852e−04 | +1.2665e−11 |
| 0.30 | 5 | vmc-ft | **−4.871661e−05** | +9.784578e−12 | −1.5852e−04 | +1.7334e−11 |
| 0.30 | 10 | ho | **−4.852497e−05** | +2.059353e−11 | −1.5790e−04 | +5.0374e−11 |
| 0.30 | 10 | vmc-ft | **−4.852497e−05** | +3.938178e−11 | −1.5790e−04 | +6.9223e−11 |

**What the correction is, exactly.** The ratio published ÷ shipped is
**3.253983147 at every one of the eight points**, and it is not a fit: it is
`1 / (LI6_B1_RANK2_TRANSFER × LI6_B1_PER_NUCLEON)` = `1 / (0.921947 × 2/6)`,
the two factors `Li6B1` applies to the deuteron table and `toy_b1` does not.
`A_zz` is **exactly linear in b₁** here because the shipped `b2` is the
Callan–Gross `2·x·b₁` — measured, `tables().b2 == 2*x*tables().b1` as doubles
at all three x — so the whole column divides by one constant. `ΔA_zz` does
**not**, because `ΔA_zz = [2 r_T − A_zz r_U] / (1 + r_U)` carries an
`A_zz`-independent term as well.

**Three readings that change, and two that do not.**

1. **The sign flip at x = 0.01 is not one row, it is all six.** The 2026-09-04
   note recorded the flip at `(0.01, 5)` `ho`. Measured here, `ΔA_zz` changes
   sign at **x = 0.01 on both C0 edges at every Q² in the table** (+8.3816e−08
   → −2.355530e−08, +3.5167e−07 → −1.769797e−07, +1.0939e−06 → −8.419210e−07
   on `ho`; the same three on `vmc-ft`). Where `2 r_T` and `−A_zz r_U` are
   comparable, shrinking `A_zz` by 3.25 tips the sum.
2. **At x = 0.10 and 0.30 the sign survives and the magnitude falls, but not
   by one factor**: ×0.318 / ×0.319 / ×0.319 at x = 0.10 (Q² = 2/5/10, `ho`)
   and ×0.405 / ×0.409 at x = 0.30 — because only part of `ΔA_zz` carries
   `A_zz`.
3. **The C0 SHAPE band on `ΔA_zz` widens rather than shrinks** at x = 0.01,
   Q² = 5: `ho` → `vmc-ft` was +3.516736e−07 → +3.229021e−07 (spread
   **2.877e−08**) and is now −1.769797e−07 → −2.100985e−07 (spread
   **3.312e−08**, +15 %) — the two edges have different `r_U`, so the
   `−A_zz r_U` term does not cancel between them.
4. **Unchanged: item 1's sign-change-in-x conclusion and item 2's
   `σ^q_U/σ^el_U` ladder**, which contain no `A_zz` at all.
5. **Unchanged: item 3's ordering.** The band still leads the whole tail at
   x = 0.01, Q² = 5 — 1.358472e−04 ÷ 1.769797e−07 = **×767.6**, the 768
   already recorded there.

**Reproduction of the ORIGINAL column, stated with its exceptions.** Every
`ho` row of §8.1c reproduces here to every printed digit, `A_zz` and `ΔA_zz`
alike, as do §8.2's five rows in full. **Three of the eight `vmc-ft` `ΔA_zz`
entries do not**: +7.9026e−08 / +3.2291e−07 / +9.7195e−07 / +2.5575e−08 as
published against **+7.902804e−08 / +3.229021e−07 / +9.720112e−07 /
+2.557553e−08** re-measured — a **2.4 × 10⁻⁵ to 6.3 × 10⁻⁵ relative**
difference in the fifth significant figure at `(0.01, 2)`, `(0.01, 5)`,
`(0.01, 10)` and `(0.10, 10)` (the other four `vmc-ft` entries reproduce).
The re-measurement is deterministic in this tree — two separate processes give
bit-identical values, and an explicit `--rc-c0-shape ho` is bit-identical to
the default — so this is a difference against the 2026-09-04 build and **its
cause is not established**. It is recorded rather than smoothed over; it is
four to five orders of magnitude below every conclusion drawn from the column,
and the `ho` edge, which carries all of them, is exact.

**Recipe** (from `env.sh`, nothing written):

```python
import lipolgen as lg
def pipe(**c0):
    cfg = lg.make_config(isotope="6Li", channel="inclusive", config=1,
                         events=2000, seed=99, rc="tensor-band", **c0)
    return lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
base = pipe(); k, s = base.dis_sampler.kernel, base.dis_sampler.s
for x, q2 in ((0.01, 2.0), (0.01, 5.0), (0.01, 10.0), (0.10, 2.0),
              (0.10, 5.0), (0.10, 10.0), (0.30, 5.0), (0.30, 10.0)):
    y = q2 / (x * s); t = k.tables(x, q2)
    azz6 = lg.azz(t.b1, t.f1, t.f2, x, y, t.b2, 0.0)              # SHIPPED 6Li b1
    azzd = lg.azz(lg.toy_b1(x, q2, t.f1, lg.B1Mode.Digitized),
                  t.f1, t.f2, x, y, None, 0.0)                    # the published column
    for p in (base, pipe(rc_c0_shape="vmc-ft")):
        m = p.rc_model
        rU = m.tail_ratio_at(x, q2, 0.0); rT = m.tail_ratio_at(x, q2, 1.0) - rU
        print(x, q2, azz6, azzd, (azz6 + 2*rT)/(1 + rU) - azz6)
```

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

> **⚠ The `A_zz` column here carries the same ×3.253983 error as §8.1c's** (deuteron b₁ instead of ⁶Li b₁), so `τ`, `w_hi − 1` and `w_lo − 1` are all high by that factor. Correct ⁶Li `A_zz` at Q² = 5: −4.528242e−04 at x = 0.01, −1.528923e−03 at x = 0.10, −4.871661e−05 at x = 0.30. Not republished here — `../run_2026-09-03/phase_B_numbers.md` §B3.2. **REPUBLISHED BESIDE IT ON 2026-09-06 — §8.2-corr immediately below; the table above stays as published.**
`w_lo + w_hi == 2.0` bit-for-bit on every event (T4a). The band peaks near
`x = 0.063` — not at the lowest x — because `δ(x)` is still rising there while
`|A_zz|` has not yet fallen: **the largest RC systematic on `A_zz` sits in the
middle of the low-x range, not at its edge.**

### 8.2-corr The same band with the SHIPPED ⁶Li b₁ — 2026-09-06

Same configuration, same `δ(x)`, same `τ` definition; `A_zz` from
`InclusiveKernel::tables` (`Li6B1(MillerB1)`) instead of the deuteron table.
**`δ(x)` is `rc_delta` alone and does not move.** The last column is the band
**half-width on `A_zz`**, `δ(x)·|A_zz|`, which is the quantity every RC-budget
comparison in this tree uses — a *different* number from `w_hi − 1` (`τ` is
≈ `A_zz`/2), and the one registry rows 6, 17 and B6 all quote.

| x | δ(x) | `A_zz` | `τ` | `w_hi − 1` | `w_lo − 1` | half-width δ·\|A_zz\| |
|---|---|---|---|---|---|---|
| 0.010 | 0.30000 | **−4.528242e−04** | **−2.264634e−04** | **−6.793901e−05** | +6.793901e−05 | 1.358472e−04 |
| 0.063 | 0.11081 | **−1.316774e−03** | **−6.588209e−04** | **−7.300143e−05** | +7.300143e−05 | **1.459067e−04** |
| 0.100 | 0.06331 | **−1.528923e−03** | **−7.650466e−04** | **−4.843711e−05** | +4.843711e−05 | 9.680016e−05 |
| 0.160 | 0.01500 | **−1.280460e−03** | **−6.406404e−04** | **−9.609605e−06** | +9.609605e−06 | 1.920691e−05 |
| 0.300 | 0.01500 | **−4.871661e−05** | **−2.435890e−05** | **−3.653835e−07** | +3.653835e−07 | 7.307492e−07 |

**The reading survives the correction, measured.** The band still **peaks at
x = 0.063**: |w_hi − 1| there is ×1.0745 its x = 0.010 value (×1.0756 on the
published column), and the half-width ×1.0740. `w_lo + w_hi == 2.0` is
unaffected — it is an identity of the construction, not of `A_zz`. Every row
of the published table above reproduces to every printed digit, so the
correction is a clean division and not a re-measurement: `A_zz` is exactly
`published ÷ 3.253983147` at all five x, while `τ` and `w_hi − 1` are only
nearly so (τ's denominator is non-linear).

**Where else these five rows are quoted, and what each site now says.**
`OPEN_ITEMS_SOLUTIONS.md` §9 carries the same table and is corrected beside its
original at the same date; the x = 0.010 half-width **1.358472e−04** is the one
`../run_2026-09-03/phase_B_numbers.md` §B3.2/§B6 and
`../run_2026-09-06/phase_A_numbers.md` §A2 already use, so this table is now
consistent with them rather than 3.25× above them.

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
| fitted `(a, α, q₀)` and χ²/ndf against the Li/Suelzle data | **NOT DONE.** `a = 1.9069 fm`, `α = 0.13822`, `q₀ = 3.0999 fm⁻¹` are the design's **unfitted starting values**, and they are shipped as such: the Suelzle–Yearian–Crannell (PR 162 (1967) 992) and Li–Sick–Whitney–Yearian (NPA 162 (1971) 583) elastic data are **not in this repository in any machine-readable form**, so there is nothing to fit against. They do reproduce their two design targets exactly (`⟨r²⟩_point = 6.0788 fm²`, first C0 zero at 3.0999 fm⁻¹, both asserted in T11). **Q1 stays open.** **PRICED 2026-09-04 (task B1), still not fitted:** the shape now runs as a two-edge band, `--rc-c0-shape ho | vmc-ft`, the second edge being the j₀ transform of the committed ANL VMC point-proton density at the *same* measured ⟨r²⟩_point — a shape a two-parameter oscillator cannot be refitted to reach. The two edges differ by ×2.5 in `F_point` at q = 2 fm⁻¹ and **flip the sign of `(1/6)σ^el_T/σ^el_U` at x = 0.10** (§8.1c). A band is a price tag, not an answer: no fit, no χ²/ndf, no data. `../run_2026-09-03/phase_B_numbers.md` §B1. |
| fitted `(q_z, b)` of `F_mag` against WS98's three `F_T` landmarks | **NOT DONE.** `q_z = 1.30 fm⁻¹`, `b = 1.85 fm`. WS98 (nucl-th/9807037) is **figures only — there is no table anywhere in the paper** (design §2.1), and a two-parameter shape cannot hit three landmarks in any case. `--rc-tail-tensor-scale 0.5 / 1 / 2` bands the NORMALISATION, not the shape. **Q10 stays open.** |
| δ(A_zz) at `x = 0.01, Q² = 5` from the band, `δ_low = 0.19` vs `0.30` | **2.7996e−04 vs 4.4204e−04** (half-widths on `A_zz`) — **the dominant RC systematic, by a factor 300 over everything below** |
| δ(A_zz) at `x = 0.01, Q² = 5` from `a_transfer_frac = 0 / 0.5 / 1` (the A = 2 → A = 6 TRANSFER price, 2026-09-04, task B6) | half-widths on `A_zz` **1.358472e−04 → 1.518818e−04 → 1.921170e−04**, i.e. the band widens by `√(1 + f²)` = **1× / 1.118× / 1.414×** (+0 % / +11.8 % / +41.4 %). Measured with **this generator's own ⁶Li b₁** (`A_zz = −4.528242e−04`), so the `f = 0` entry is the **corrected** 1.3585e−04 of §B3.2 and **not** the 4.4204e−04 of the row above. Against the whole tail's `ΔA_zz = −1.769797e−07` at the same point the band is **×768 / ×858 / ×1086** — the factor-1300 headline elsewhere in this file used the same ×3.253983-too-large `A_zz`; the ordering is unchanged, the factor is **768**, not 1300. Default **0.0**, and the shipped npz is byte-identical to the pre-B6 one (51 columns + `meta`, both `--rc tensor-band` and `--rc off`). `../run_2026-09-03/phase_B_numbers.md` §B6 |
| δ(A_zz) at `x = 0.01, Q² = 5` from `fq_scale = 0` vs `2` | `ΔA_zz` **+7.54154e−07 → −4.9925e−08**, spread **8.041e−07** |
| δ(A_zz) at `x = 0.01, Q² = 5` from `tail_tensor_scale = 0.5` vs `2` | `ΔA_zz` **+3.59193e−07 → +3.23118e−07**, spread **3.61e−08** |
| δ(A_zz) at `x = 0.01, Q² = 5` from `qe_suppression = 0` vs `1` | `ΔA_zz` **+1.83644e−07 → +3.51674e−07**, spread **1.680e−07** |
| δ(A_zz) at `x = 0.01, Q² = 5` from `qe_kf_gev = 0` (S = 1) vs `0.169` | `ΔA_zz` **+5.39491e−07 → +3.51674e−07**, spread **1.878e−07** — the Pauli suppression is a *larger* knob than `qe_suppression = 0 → 1` at this x, which is why it is the default rather than a band edge |

### 8.3-corr The six `A_zz`-carrying rows of that table, recomputed — 2026-09-06

The five `ΔA_zz` rows and the `δ_low` row above all carry the ×3.253983
deuteron-b₁ `A_zz` (§8.1c-corr). They are **left as published** and recomputed
here with the shipped ⁶Li b₁, same point (x = 0.01, Q² = 5), same knob values,
same `ho` C0 edge unless the row names another.

| row | as published | **with the shipped ⁶Li b₁** | spread, published → corrected |
|---|---|---|---|
| band, `δ_low = 0.19` vs `0.30` (half-widths on `A_zz`) | 2.7996e−04 vs 4.4204e−04 | **8.603659e−05 vs 1.358472e−04** | 1.6208e−04 → **4.981061e−05** |
| `fq_scale` 0 vs 2 | +7.54154e−07 → −4.9925e−08 | **+2.255019e−07 → −5.785810e−07** | 8.041e−07 → **8.041e−07** (unchanged) |
| `tail_tensor_scale` 0.5 vs 2 | +3.59193e−07 → +3.23118e−07 | **−1.694451e−07 → −2.055947e−07** | 3.61e−08 → **3.615e−08** |
| `qe_suppression` 0 vs 1 | +1.83644e−07 → +3.51674e−07 | **−2.286503e−07 → −1.769797e−07** | 1.680e−07 → **5.167e−08** |
| `qe_kf_gev` 0 (S = 1) vs 0.169 | +5.39491e−07 → +3.51674e−07 | **−1.192240e−07 → −1.769797e−07** | 1.878e−07 → **5.776e−08** |
| `c0_shape` `ho` vs `vmc-ft` (§9's row) | +3.51683e−07 → +3.22911e−07 *(§9 as printed; this tree re-measures the same two settings at +3.516736e−07 → +3.229021e−07 — see the exception note)* | **−1.769797e−07 → −2.100985e−07** | 2.877e−08 → **3.312e−08** |

**Five of the six "as published" pairs reproduce exactly** from this tree with
the deuteron table in place of the kernel's b₁, to every digit they print —
`fq_scale`, `tail_tensor_scale`, `qe_suppression`, `qe_kf_gev` and both
half-widths of the `δ_low` row (0.19 × 1.473482e−03 and 0.30 × 1.473482e−03).
**The `c0_shape` pair does not**: §9 prints +3.51683e−07 → +3.22911e−07 where
the same two settings here give +3.516736e−07 → +3.229021e−07 — 2.7 × 10⁻⁵ and
2.9 × 10⁻⁵ relative away, in the fifth significant figure — while its
**spread**, 2.8772e−08, reproduces exactly. It is the same `vmc-ft` exception
§8.1c-corr records, and its cause is likewise not established; the `ho` value
in that published pair also disagrees with the four other rows' own `ho` value
inside §9, which is what makes it a property of that measurement rather than
of the edge.

**Which spreads move, and why — the rule is one line.** A knob that leaves
`r_U` alone has an `A_zz`-independent spread (`Δ(ΔA_zz) = 2 Δr_T/(1 + r_U)`);
a knob that moves `r_U` does not, because the `−A_zz r_U` term then differs
between its two settings. So `fq_scale` (which touches only `F_q`, i.e. `r_T`)
is **unchanged to five digits** (8.040793e−07 vs 8.040830e−07), `tail_tensor_scale` moves by +0.21 % — the third significant figure (3.607508e−08 → 3.614951e−08),
and the three knobs that scale the quasi-elastic tail or the elastic C0 shape —
`qe_suppression`, `qe_kf_gev`, `c0_shape` — move by ×0.31 / ×0.31 / ×1.15.

**What this does to the ORDER of the RC budget, which is the reading the table
exists for.** Within §8.3's own four `ΔA_zz` knobs the order is **unchanged**:
`fq_scale` ≫ `qe_kf_gev` > `qe_suppression` > `tail_tensor_scale`. It changes
once §9's two extra rows are included, and it changes for a bookkeeping reason
worth naming: §9 already carried `qe_tensor_scale` **corrected** (1.162e−07,
measured on the ⁶Li b₁ in 2026-09-04's task B3) beside four rows that were
not, so its published ladder compared one corrected number against four
uncorrected ones. All six on the same footing:

**`fq_scale` 8.041e−07 ≫ `qe_tensor_scale` 1.162e−07 > `qe_kf_gev` 5.776e−08 >
`qe_suppression` 5.167e−08 > `tail_tensor_scale` 3.615e−08 > `c0_shape`
3.312e−08** — `qe_tensor_scale` moves from **fourth to second**, not because it
grew but because the three below it shrank by ×0.31. `c0_shape` stays the
smallest of the six, which is what §9's own sentence claims.

**And the band still dominates all of them** — the ordering the whole section
turns on. The published row's *"by a factor 300 over everything below"* does
not reproduce as either well-defined ratio even on its own column (the
half-width over the largest knob spread is **549.7**; the 0.19 ↔ 0.30 spread
over it is **201.6**), so it is restated here as the two measured numbers:
with the shipped ⁶Li b₁ the half-width is **168.9×** the largest knob spread
and the 0.19 ↔ 0.30 spread is **61.9×** it. **The band is still the dominant
RC systematic on both readings**, by two orders of magnitude and one order
respectively.

**The `fq_scale` band CHANGES THE SIGN of `ΔA_zz`** (+3.9e−06 at 0, −9.6e−07 at
2 — *those two numbers are older than the table above them: they reproduce
under neither b₁ choice, and they predate the 2026-09-03 per-nucleon fix.
Flagged 2026-09-06, left in place; the reproducing pairs are the `fq_scale`
row of §8.3 and of §8.3-corr, and the SIGN CHANGE — the claim of this
sentence — holds on both of them*) and is **not** symmetric about the nominal — `σ^el_T` is *quadratic* in
`F_q`. That is exactly why the design forbids rescaling one run and why T12 fits
a quadratic through three runs instead of asserting linearity.

**Clipping.** **Nothing clips any more at the default `tail_model = TPeak`.**
Before the 2026-09-03 per-nucleon fix
21.95 % of the `y > 0.9` nodes exceeded `tail_max = 10` and one event in 2000
hit the ceiling in the default CLI run; dividing the elastic tail by the missing
`A = 6` removed all of it. On the opt-in `--rc-tail-model t-peak+ll` the
**nodes** clip again — **0.36 %** globally and **4.62 %** of the `y > 0.9`
band on this 101 × 77 grid (1.01 % / 5.85 % on T8(d)'s coarser 40 × 24 test
grid) — while the **event** count stays **0** at 2 000 and at 200 000 events,
because no accepted cell centre lands in the nodes that clip. The two
statistics diverging is precisely why they are quoted separately. The `y > 0.9` corner is still where `Y₊ ~ 1/(1−y)`
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
* The **polarised** quasi-elastic tail is **not computed anywhere**, and
  since 2026-09-04 it is *priced* by the opt-in stand-in
  `--rc-qe-tensor-scale` (default 0.0 = this document's numbers, unchanged),
  which lends it the elastic tail's own tensor fraction — a borrowed
  magnitude, not a derived bound (`../run_2026-09-03/phase_B_numbers.md` §B3).
  It is neglected by default citing Zhou et al., PRL 82 (1999) 687 -- a
  *deuteron* statement -- and at `x = 0.30` the tail is **99.9 %**
  quasi-elastic, which is what makes that zero the largest unpriced piece.
* **`rc_tail` is the t-PEAK ONLY and is therefore a LOWER BOUND.** Measured
  against the leading-log s-/p-peaks (T8(c), T8(d)(i)): at `Q² ≥ 20 GeV²` the
  "> 99 %" this bullet carried until 2026-09-04 is **EVENT-WEIGHTED ONLY** —
  the window mean moves **+0.61 %** — and is **false cell by cell**:
  over `Q² ≥ 20 GeV²` with no `y` cut, **346 of 1399 accepted cells (24.7 %,
  29.6 % of that window's cross section) differ by more than 1 %**, **318 of
  them at `y < 0.1`**, worst **×6444**. The per-cell statement that survives is
  `Q² ≥ 20 GeV²` **and** `0.15 ≤ y ≤ 0.7` (660 cells, worst **0.55 %**); §8.1
  above carries the same correction and the mechanism. The lower bound
  elsewhere is unchanged: **low by a factor 4.4 at the HERMES deuteron point**
  and by 3.3× on the quasi-elastic piece at the `Q² ≈ 4 GeV²` corner of the
  generator window.
* The **tagged band is clamped** (`RcOptions::band_tau_max = 1`). `tau_tag`
  itself is unbounded at the nodes of `n_M(k, c)`; the clamp is a *choice*, the
  clipped fraction is reported, and `n_M → 0` is exactly where the
  fractional-rescale ansatz stops meaning anything.
