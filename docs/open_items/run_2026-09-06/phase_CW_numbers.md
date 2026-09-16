<!-- SPDX-License-Identifier: GPL-3.0-or-later -->
# Run 2026-09-06 — Phase CW measurements: the tagged S–D sign fix

**This phase is a BUG FIX, the one exception to "no default moves" in this
run.** `docs/benchmarking/07_cw_sign_investigation.md` returned the tagged
sector's S–D interference sign **refuted at blocker severity**, with three
independent derivations and the verdict *certain*: `TaggedModel::build_amp2`
summed the partial waves as ψ_L where a partial-wave amplitude needs
**φ_L = i^L ψ_L**. Parity fixes L mod 2 inside one channel, so the common
`i^(L mod 2)` is an unobservable global phase and the observable **relative**
phase is real — `(-1)^floor(L/2)`, i.e. **+1 on L = 0, +1 on L = 1, −1 on
L = 2** — and it was missing entirely.

Why it is a bug and not a registry decision, restated from the investigation
§7.3 and re-verified here: it disagrees with a **published formula**
(Cosyn–Weiss II Eq. 6.12, by up to **2.740** in an asymmetry whose whole range
is `[−2, 1]`), with the repository's **own deuteron quadrupole sign gate**
(the shipped amplitude gives Q_d = −0.3036 fm² against the file's measured
+0.2697), and with the repository's **own b₁ sector**, which applies the same
phase explicitly in `ClusterPartialWave::from_vmc` / `::from_uw` and documents
it as LOAD-BEARING.

**The fix is one line** (`src/core/tagged.cpp` `build_amp2`), plus the two doc
comments that carried the false inference, four re-pinned tests, a fifth
re-pin the investigation did not enumerate (T25's dilutions, §7 below), and a
**deliberate re-pin of `validation/reference/tagged.json`** (§8). The stored
`ψ₂ = +W` convention is untouched: `rad_[2]` and `Wave::vmc->psi()` are not
negated, and the two tests that correctly assert that convention still pass.

---

## 0. Suites and gates

| | before | after |
|---|---|---|
| C++ `build/lipolgen_tests` | 404 cases, **404 passed / 0 failed / 1 skipped**, 17 240 384 assertions | 404 cases, **404 passed / 0 failed / 1 skipped**, **17 240 391** assertions |
| Python `pytest python/tests -q` | **1004 passed / 127 skipped** | **1004 passed / 127 skipped** |
| docs gate (strict) | 1220 refs / **95** ranges / 7 external / **0 broken** / 6 allow-listed; +19/115 `SPIN32_FINITE_GAMMA.md`; +8 external `PYTHIA_BRIDGE.md` | 1220 refs / **96** ranges / 7 external / **0 broken** / 6 allow-listed; same sidecars |
| SPDX | **101/101** | **102/102** |
| `validation/reference/` | — | **only `tagged.json` and `_manifest.json` moved**; the other seven are byte-identical (SHA-256 verified, §8.3) |

The +7 assertions are the rewritten Cosyn–Weiss gate's net. The +1 range is
one new citation added to `PHYSICS_CHANNELS.md`'s A_zz^tag row
(`docs/open_items/vmc_reconciliation.md:180-188`, which now carries a
correction box). The +1 SPDX file is
`validation/repin_tagged_from_lipolgen.py`. Eight moved doc-block references
and one moved-symbol reference were re-anchored with `--fix`; three dropped
and one moved range fingerprint with `--record-ranges`. **Not one cited block
was edited to make the gate pass** — every re-anchor is a line shift caused by
the comment and correction-box insertions.

Measurement method: a scratch probe compiled read-only against the shipped
`build/libLiPolGenCore.so`, run once against a copy of the **pre-fix** library
and once against the **post-fix** one, same binary source, same data
directory. Every "before" number below reproduces the investigation's own
"shipped" column and every "after" number reproduces its "i^L-applied"
column, which is what makes this a like-for-like measurement of one phase.

**Second instrument, and it needs no second library.** Because `build_amp2`
now applies `(-1)^floor(L/2)`, **negating the stored L = 2 table on the fixed
library reproduces the pre-fix amplitude exactly** — the two factors cancel.
So a "before" column can be measured in-process, on one binary, by installing
`−radial_table(2)` as the channel's L = 2 `VmcRadial`. Its validation is that
it reproduces this tree's **own published pre-fix numbers cell for cell**:
all fifteen cells of `vmc_reconciliation.md:180-188` (§2), the four channels'
dilutions and `p2_moment`s (§4, §7.2), the CW identity residual
**2.740499e+00** (§6) and the three TABLE II rows **+0.617024 / −0.318307 /
−1.656113**. Every number in the 2026-09-06 documentation pass was measured
with it, independently of the probe above, and the two agree.

**Gate deltas of the documentation pass** (the second pass of §10, which
touched no `src/` and no `include/`): the docs gate went **0 → 63 broken →
0**. Sixty were plain line citations that `--fix` re-anchored, all of them
shifted by the two source insertions this pass made — a help clause in
`python/lipolgen/cli.py` and a docstring paragraph in
`python/lipolgen/configs.py`. The remaining three `--fix` could not resolve;
they were re-pointed **by content**, checked byte-for-byte against their
recorded fingerprints, and then re-recorded — `docs/USAGE.md:2530-2537 → :2534-2541`,
`python/lipolgen/cli.py:1346 → :1354` (`rc_model`) and
`python/lipolgen/cli.py:1109 → :1117` (`set_pythia_hadronizer`). **The
fingerprint file proves no cited content changed**: comparing it against
`HEAD`, the number of entries whose recorded `sha256`/`first`/`last` differs
is **zero** — all fifteen moves are pure re-anchoring. Suites and SPDX are
unchanged by this pass (404 / 17 240 391 / 1 skipped; 1004 / 127; 102/102).

---

## 1. The line

`src/core/tagged.cpp`, inside `build_amp2`'s per-wave angular prefactor loop:

```cpp
    for (std::size_t w = 0; w < rad.size(); ++w) {
      // phi_L = i^L psi_L.  Parity fixes L mod 2 inside one channel, so the
      // common i^(L mod 2) is an unobservable GLOBAL phase and the observable
      // RELATIVE phase is real: (-1)^floor(L/2).  Identical rule to
      // ClusterPartialWave::from_vmc / from_uw (src/core/b1_nuclear.cpp:249,
      // :270), which is why the b1 and tagged sectors then agree.
      const double phase = ((l_of[w] / 2) % 2 == 0) ? 1.0 : -1.0;
      for (std::size_t ic = 0; ic < nc; ++ic) {
        ang[w][ic] = phase * cg_of[w] * theta_lm(l_of[w], ml_int, c_[ic]);
      }
    }
```

**Verified on the grid, not asserted:** `L = 0 → +1`; `L = 1 → +1`, and ⁷Li
α-tag is measured **bit-identical** across the fix (§3); `L = 2 → −1`, and
every spin-1 channel's S–D interference flips (§2, §4, §5).

The two doc comments that carried the false inference —
`include/lipolgen/tagged.hpp` and `src/core/tagged.cpp`, the same sentence in
both, *"CDKS fix φ_L = i^L ψ_L with φ₂ = −W and U, W ≥ 0 at low k, so
ψ₂ = +W"* — now say what is true of each object: **ψ₂ = +W is the STORED
convention and is right; `build_amp2` consumes φ, not ψ, and applies
`(-1)^floor(L/2)` itself.**

---

## 2. The headline — `A_zz^tag(k)`, ⁶Li, acceptance-weighted

`azz_tensor_curve_weighted` at `default_configs("6Li")[1]` (10 × 99.5 GeV/u),
`yr_optics("6Li", ..., high_acceptance = true)`, `n_phi = 32`. Cells at the
grid points nearest the quoted k.

| k [GeV] | Hulthén β=0.30 **before** | **after** | Hulthén P_D=VMC **before** | **after** | VMC AV18 **before** | **after** |
|---|---|---|---|---|---|---|
| 0.1979 | +0.8450 | **−1.2069** | +0.5108 | **−0.6001** | +0.4518 | **−0.5191** |
| 0.2495 | +0.6551 | **−1.0639** | +0.4639 | **−0.5767** | +0.2604 | **−0.2899** |
| 0.3012 | +0.5226 | **−0.9533** | +0.4296 | **−0.5615** | +0.1819 | **−0.1993** |
| 0.4001 | +0.3325 | **−0.7267** | +0.3435 | **−0.4839** | +0.0779 | **−0.0822** |
| 0.4990 | +0.2349 | **−0.5851** | +0.2840 | **−0.4214** | −0.1302 | **+0.1170** |

The "before" columns are `docs/open_items/vmc_reconciliation.md:180-188`
exactly, and **the "after" VMC column is that document's `vmc-flipD` control
column, digit for digit** (−0.5191, −0.2899, −0.1993, −0.0822, +0.1170). The
tree's designated *control on the S–D relative sign* was the physical answer.
The old physics column is now the control: negating the D table on the fixed
code gives back +0.4518 / +0.2604 / +0.1819 / +0.0779 / −0.1302.

---

## 3. ⁷Li α-tag — must not move, and does not

One L = 1 wave: nothing to interfere with, so the common `i` is a global phase
and every observable is invariant. Measured at 17 significant digits, Hulthén
and VMC, both before and after:

| quantity | before | after |
|---|---|---|
| ⁷Li-H `n_of_kc(+3/2)` Σ over the grid | `6.74976676548482646e+04` | `6.74976676548482646e+04` |
| ⁷Li-H `n_of_kc(+3/2)` cell [100] | `6.16484018725258487e-02` | `6.16484018725258487e-02` |
| ⁷Li-H `n_of_kc(+3/2)` cell [5000] | `2.22158508416940892e+00` | `2.22158508416940892e+00` |
| ⁷Li-H `n_of_kc(+1/2)` Σ | `6.74903440865036246e+04` | `6.74903440865036246e+04` |
| ⁷Li-V `n_of_kc(+3/2)` Σ | `1.20278741037640575e+05` | `1.20278741037640575e+05` |
| ⁷Li-V `n_of_kc(+1/2)` Σ | `1.20265690666996117e+05` | `1.20265690666996117e+05` |
| `p2_moment(±3/2)` | −0.1999349117 | −0.1999349117 |
| `p2_moment(±1/2)` | +0.1998480944 | +0.1998480944 |
| `p2_moment_mixture_uniform` | −4.3408664113e−05 | −4.3408664113e−05 |
| `population_integrated(M)`, all four M | (1, 0) / (0.6666304957, 0.3333695043) / … | identical |
| `norm(±3/2)` / `norm(±1/2)` | 1.0001232231 / 1.0000147087 | identical |
| `vector_dilution` | 1.0000000000 | 1.0000000000 |

**Bit-identical, every entry.** ⁷Li's block of `tagged.json` was therefore
left as polligen wrote it (§8.2).

---

## 4. `p2_moment` per M, every spin-1 channel

The angle-integrated ⟨P₂(cos θ_k)⟩ per ion projection. Every one flips sign;
the magnitudes grow because the corrected S–D term adds to the D² term
instead of cancelling against it.

| channel | M | before | after | ratio |
|---|---|---|---|---|
| ⁶Li α-tag Hulthén | +1 / −1 | +0.0448133305 | **−0.0622510957** | ×1.389 |
| ⁶Li α-tag Hulthén | 0 | −0.0897861433 | **+0.1243457311** | ×1.385 |
| ⁶Li α-tag VMC AV18 | +1 / −1 | +0.0038621121 | **−0.0078394851** | ×2.030 |
| ⁶Li α-tag VMC AV18 | 0 | −0.0078869596 | **+0.0155163085** | ×1.967 |
| deuteron control Hulthén | +1 / −1 | +0.0320226969 | **−0.0411257151** | ×1.284 |
| deuteron control Hulthén | 0 | −0.0642064787 | **+0.0820914169** | ×1.279 |
| deuteron control AV18 | +1 / −1 | +0.0198147059 | **−0.0314365982** | ×1.587 |
| deuteron control AV18 | 0 | −0.0397915293 | **+0.0627120399** | ×1.576 |
| ⁷Li α-tag (both wave sources) | all four | (see §3) | **unchanged to the bit** | ×1 |

`p2_moment(M=0) = −2 p2_moment(M=±1)` holds exactly before and after — it is
CG algebra, not a sign statement.

---

## 5. `A_zz^wf(k)` at the θ_k ≈ 0 and θ_k ≈ 90° cells

The unweighted cell curves, at |c| = 0.989583 and |c| = 0.010417.

**Deuteron control, Hulthén β = 0.30:**

| k [GeV] | θ≈0 before | θ≈0 after | θ≈90° before | θ≈90° after |
|---|---|---|---|---|
| 0.0990 | +0.344326 | **−0.378651** | −0.177629 | **+0.195336** |
| 0.1979 | +0.827412 | **−1.111630** | −0.426841 | **+0.573462** |
| 0.2968 | +0.967648 | **−1.575777** | −0.499185 | **+0.812903** |
| 0.4001 | +0.937986 | **−1.778248** | −0.483883 | **+0.917353** |
| 0.5979 | +0.837952 | **−1.891958** | −0.432278 | **+0.976013** |
| 0.9979 | +0.743817 | **−1.928114** | −0.383716 | **+0.994664** |

**Deuteron control, AV18 (`ClusterWaveSource::VmcAV18`) — CW's own wave
function.** The last column is Cosyn–Weiss II Eq. (6.12) evaluated
independently at the same cell:

| k [GeV] | f₂/f₀ | θ≈0 before | θ≈0 after | CW Eq. (6.12) |
|---|---|---|---|---|
| 0.0990 | +0.1142 | +0.296342 | **−0.321268** | −0.321268 |
| 0.1979 | +0.4654 | +0.875862 | **−1.220872** | −1.220872 |
| 0.2968 | +1.3942 | +0.658141 | **−1.937694** | −1.937694 |
| 0.4001 | (past the U node) | −0.799965 | **−1.130776** | −1.130776 |
| 0.5979 | −1.6947 | −1.918141 | **+0.480791** | +0.480791 |
| 0.9979 | −0.7426 | −1.656113 | **+0.967340** | +0.967340 |

**⁶Li α-tag, Hulthén β = 0.30:** θ≈0 goes +0.406856 / +0.901583 / +0.955676 /
+0.859119 / +0.707209 / +0.589062 → **−0.456153 / −1.287665 / −1.727568 /
−1.877767 / −1.934249 / −1.935204** at the same six k.

**⁶Li α-tag, VMC AV18:** θ≈0 goes −0.443879 / +0.481984 / +0.336227 /
+0.201334 / −1.910663 → **+0.397127 / −0.553864 / −0.368844 / −0.212442 /
+0.449030** at the first five k. (The k = 0.9979 cell is NaN on this channel
before and after — the VMC tables end below it, and the acceptance-weighted
curve of §2 is the quantity that is actually quoted there.)

---

## 6. Cosyn–Weiss TABLE II on the AV18 control, with the corrected code

CW's TABLE II is quoted *"based on the AV18 radial wave functions"*, so the
rows run on `deuteron_channel(..., ClusterWaveSource::VmcAV18)` — the
like-for-likeness failure `06` flagged, and the reason the gate moved off
Hulthén.

| row | k, θ_k | CW `A_T∥` | **measured after** | before |
|---|---|---|---|---|
| 1 | 0.30 GeV, 0 | **−2** | **−1.937124** | +0.617024 |
| 2 | 0.30 GeV, π/2 | **+1** | **+0.999313** | −0.318307 |
| 3 | 1 GeV, 0 | **+1** | **+0.967340** | −1.656113 |

The residuals from CW's exact −2 / +1 / +1 are the finite cell centres:
`c = ±0.989583` rather than ±1 (P₂ = 0.968913, not 1) and `k = 0.30115` rather
than the exact `f₂/f₀ = √2` point. Through the P₂ factorization, which is
exact to **1.2 × 10⁻¹⁴**, the extremes come out at CW's values:

| quantity | measured after | CW | before |
|---|---|---|---|
| `A_zz^wf/P₂` envelope minimum | **−1.999864** at k = 0.296849 (f₂/f₀ = +1.3942 → √2) | −2 | −1.999999 at k = 0.649508 |
| `A_zz^wf/P₂` envelope maximum | **+0.999997** at k = 1.036573 (f₂/f₀ = −0.7056 → −1/√2) | +1 | +0.999991 at k = 0.236639 |
| `A_zz^wf/P₂` at k = 0.30, band midpoint | **−1.999276**, spread 1.2e−14 | −2 | +0.636821 |
| envelope min × (−1/2) | **+0.999932** | +1 (the θ = π/2 row) | — |
| k where f₂/f₀ = +√2 | **0.298121 GeV** | "0.30 GeV" | (unchanged: a wave-function property, not an amplitude one) |
| k where f₂/f₀ = −1/√2 | **1.034872 GeV** | "≈ 1 GeV" | (unchanged) |
| A_zz over the whole AV18 grid | **[−1.937694, +0.999606]** | "[−2, 1]" | [−1.937824, +0.999674] |

**And the identity, over every one of the 280 × 96 cells of the grid:**

    max | A_zz^wf(k,c) − [(2f₀ + f₂/√2)(f₂/√2)/(f₀²+f₂²)] (1 − 3cos²θ_k) |

| channel | before | after |
|---|---|---|
| deuteron AV18 | **2.740499e+00** | **8.881784e−16** |
| deuteron Hulthén | **2.740499e+00** | **8.881784e−16** |

That is the whole of CW Eq. (6.12), to machine precision, on both wave
functions — and it settles the mapping question at the same time: the
coefficient is **`A_T∥ = +1 · A_zz^wf`**, not −2. CW's "−2" is the value of
`A_zz^wf` at θ_k = 0 where the Λ = ±1 densities have a node (their Eq. 6.13),
not a conversion factor; the old gate applied it a second time.

---

## 7. What does not move

### 7.1 The ⁶Li α-tag rate — unchanged, exactly

The spin-blind rate is the M-average of `n_M`, and `Σ_M n_M` is isotropic by
CG completeness: the S–D cross term cancels in the sum identically. Measured
as the acceptance-weighted, k²-weighted accepted fraction at the same optics
as §2:

| channel | before | after |
|---|---|---|
| ⁶Li α-tag Hulthén β=0.30 | `0.024675932148827` | `0.024675932148827` |
| ⁶Li α-tag Hulthén P_D=VMC | `0.017299520422530` | `0.017299520422530` |
| ⁶Li α-tag VMC AV18 | `0.033810227625842` | `0.033810227625842` |

**Unmoved to ≤ 4.0e−16 relative — 1–2 ulp of double summation — and the VMC
channel is bit-identical.** The unnormalized accepted rates are likewise
unmoved (`4.383238481169e+01` Hulthén, `3.072947309806e+01` Hulthén at the VMC
P_D, `6.005768857484e+01` VMC; residuals ≤ 2.3e−16). *(The VMC fraction was
first published here as `0.033810227625845`; independently re-measured
2026-09-06 it is `0.033810227625842`, 3 ulp lower, and that value is the one
that reproduces the unnormalized `6.005768857484e+01` above exactly. Nothing
physical turns on it — quote no more than `0.0338102276258`, since digits past
the 13th depend on the summation order.)* What changed at
event level is only the **correlation** the sampler draws — `sample_kc` draws
from `|A_{m_S}(M; k, c)|² k²`, so the joint (M, m_S, cos θ_k) correlation
flips while the rates do not. That is the observable.

Consequently `USAGE.md`'s tag-fraction rows (0.0249 / 0.0348 YR
high-acceptance, 0.2530 / 0.2485 tagging optics, 0.9730 / 0.9981 for ⁷Li),
⟨k⟩, P(k > ·) and P_D are all unchanged and were left alone.

> **CORRECTION, 2026-09-06 (verification pass) — this paragraph is false for
> the ⁶Li rows, and they have since been REGENERATED.** ⟨k⟩, P(k > ·), P_D and
> every ⁷Li row are unchanged (⁷Li is bit-identical). The ⁶Li tag-fraction rows
> are 40 000-event SAMPLES at the `tensor-thirds` fill, and while the
> EXPECTATION of that category average is invariant — the uniform-M mix
> accepted fraction is unmoved to all 15 digits, §7.1's own table, and the
> equal-thirds average reproduces it to **+3.1e−6 (a per-M norm-residual effect — Σ n_M k² differs by 5e−5 between M = 0 and ±1 — not summation order; the pooled category average equals the uniform mix to 2.2e−16 on both builds)** — the SAMPLE is re-drawn,
> because ~77 % of events take a different (k, cos θ_k). Re-running the rows'
> own command on the fixed library
> (`validation/vmc_tag_fractions.py --events 40000 --configs 0,1,2`, seed
> 20260829) prints **0.0301 / 0.0264 / 0.0279** (⁶Li YR high-acceptance,
> configs 0/1/2) and **0.3451 / 0.2551 / 0.3145** (⁶Li tagging) against
> 0.0286 / 0.0249 / 0.0266 and 0.3410 / 0.2530 / 0.3115 before — at config 1
> that is 1.4σ (YR) and 0.7σ (tagging) of the difference of two independent
> 40 000-event samples, binomial σ_diff 1.1e−3 and 3.1e−3 — and every ⁷Li row
> is identical digit for digit. `USAGE.md`'s table and
> `docs/open_items/vmc_reconciliation.md` are regenerated on the fixed build
> and say so; leaving them would have left published digits their own
> documented command no longer reproduces.
>
> **And a tag fraction at a single tensor-polarised fill is not invariant at
> all** — it is the same acceptance-weighted integral as A_zz^tag. Model
> predictions at the same optics, pre-fix → post-fix: ⁶Li Hulthén
> `tensor-thirds` categories **0.028127 → 0.018403** (azz±) and
> **0.017775 → 0.037222** (azz0); the `--pz 0.7` max-entropy ladder of
> `helicity-flip` **0.027030 → 0.020396** (−24.5 %), ⁶Li VMC 0.035296 →
> 0.032155, deuteron Hulthén 0.025802 → 0.021492, deuteron AV18 0.022175 →
> 0.019407. The CLI's `--cluster-wave` SIGN NOTE, which said "tag fractions …
> are unaffected", now carries this.

### 7.2 The tagged tensor dilution — unchanged as physics, 2.5e−6 as arithmetic

Angle-integrated, so L-orthogonality kills the cross term. What moves is the
**96-cell midpoint quadrature residual** of `∫Θ₀Θ₂ dc = 0`, whose sign the
phase flips — nothing physical.

| channel | quantity | before | after | relative |
|---|---|---|---|---|
| ⁶Li α-tag Hulthén | `tensor_dilution` | 0.9219467089 | 0.9219489770 | **+2.5e−06** |
| ⁶Li α-tag Hulthén | `vector_dilution` | 0.8699393995 | 0.8699431789 | +4.3e−06 |
| ⁶Li α-tag VMC AV18 | `tensor_dilution` | 0.9825758170 | 0.9825758723 | +5.6e−08 |
| ⁶Li α-tag VMC AV18 | `vector_dilution` | 0.9709659943 | 0.9709660865 | +9.5e−08 |
| deuteron Hulthén | `tensor_dilution` | 0.9594880739 | 0.9594888782 | +8.4e−07 |
| deuteron Hulthén | `vector_dilution` | 0.9324947691 | 0.9324961093 | +1.4e−06 |
| deuteron AV18 | `tensor_dilution` | 0.9481456181 | 0.9481463394 | +7.6e−07 |
| deuteron AV18 | `vector_dilution` | 0.9135947766 | 0.9135959786 | +1.3e−06 |
| ⁷Li α-tag | `vector_dilution` | 1.0000000000 | 1.0000000000 | 0 |

**`LI6_B1_RANK2_TRANSFER` = 0.921947 is unaffected** — it is pinned to the ⁶Li
Hulthén tensor dilution at 1e−4 and the move is 2.5e−6. So are the effective
polarizations, the tag fractions and the whole inclusive sector.

> **CORRECTION, 2026-09-06 (verification pass): strike "the tag fractions" from
> that sentence.** `LI6_B1_RANK2_TRANSFER`, the effective polarizations and the
> inclusive sector are unaffected as stated (the constant's quadrature reading
> moves 0.9219467089 → 0.9219489770, +2.46e−6 against a 1e−4 pin). A tag
> fraction, however, is invariant only for a spin-blind or category-averaged
> fill: it is an acceptance-weighted integral over n_M, so a tensor-polarised
> fill moves it with A_zz^tag — up to −24.5 % at the CLI's own default fill.
> Measured numbers in the box under §7.1.

**One consequence the investigation did not enumerate.** `T25 the deuteron
control channel, Hulthen against AV18` (`tests/test_tagged.cpp`) pins these
four dilutions at **rtol 1e−6**, and the two vector ones move by 1.4e−6 and
1.3e−6 — just past it. All four are re-pinned to the measured post-fix values
at the same tolerance, with the pre-fix values kept in the comment beside
them: vector **0.932494769 → 0.932496109** and **0.913594777 → 0.913595979**,
tensor **0.959488074 → 0.959488878** and **0.948145618 → 0.948146339**. The
two `v/h` RATIOS the test also pins (−0.020268, −0.011822 at rtol 1e−3) are
unchanged to every digit they carry.

### 7.3 Everything outside `tagged`

**No b₁ number moves.** `b1_nuclear` builds its own amplitudes through
`ClusterPartialWave::from_vmc`, which already applied `φ_L = i^L ψ_L`, and
only *reads* the tagged channel's `Wave::vmc` tables, which the fix leaves
untouched. Verified: `validation/reference/b1_default_li6.json` is
**byte-identical** (SHA-256 `d7bd8ce6…`), the CDKS shape gate,
`alpha_d_quadrupole_fm2` and `LI6_B1_RANK2_TRANSFER` all pass unchanged. So do
the inclusive sector, the coherent channel, the RC sector, `boost_spectator`
and every `spectator.json` / `coherent.json` / `xsec.json` / `spin.json` /
`beams.json` / `bookkeeping.json` entry (§8.3).

**The two tests that correctly assert the stored ψ convention still pass, by
construction:** `"tagged: the VMC S-D interference flips sign below the S
node"` (checks `radial_table(2)` signs about the S node) and T25's
`waves[1].vmc->psi()[1] > 0`. The fix is in the amplitude sum, not the radial
tables — which is exactly why it was done there.

---

## 8. The reference re-pin — deliberate, documented, and only `tagged.json`

### 8.1 polligen cannot be the source any more

`validation/README.md` says `tagged.json` is regenerated with
`validation/dump_polligen_reference.py`, which dumps the Python generator
`polligen`. **`polligen.tagged._amp2_table` (tagged.py:243-248) carries the
same missing `i^L`** — it sums `self._rad[w.l] * cg * theta_lm(...)` with no
phase, exactly the pre-fix C++. Regenerating from it would re-bake the refuted
sign and the rtol-1e-12 gate would go on certifying the bug.

So the two **spin-1** model blocks are re-pinned from the **fixed C++
library**, through the installed pybind11 module — the same code the tests
call, the pattern `validation/dump_b1_default_li6.py` already established. New
script: **`validation/repin_tagged_from_lipolgen.py`**, which states all of the
above in its docstring. The file now carries its own provenance keys:

    "provenance": "LiPolGen post-fix, formerly polligen"
    "provenance_note": "The li6_alpha and deuteron `model` blocks are dumped
      from the FIXED LiPolGen C++ library ... polligen's tagged._amp2_table
      omits the i^L partial-wave phase and therefore carries the INVERTED S-D
      interference sign ..."

`validation/README.md` records it in full (it is generated from
`dump_polligen_reference.py`'s own `doc()` text, so the record survives a
regeneration), and `_manifest.json`'s `generators` map now names the re-pin
script as `tagged.json`'s owner. `dump_polligen_reference.py` itself was
taught to **carry those two blocks through** rather than overwrite them, and
to refuse loudly if the file is missing — so running it can no longer
reintroduce the bug. The two scripts are idempotent against each other:
running either then the other leaves `tagged.json` byte-identical (verified,
worst move `0.000e+00`).

One binding was added to make this path possible:
`TaggedModel.struck_populations` in `python/bindings.cpp` (read-only,
`(n_mS, nk, nc)`, mirroring `n_of_kc_table`). Without it the dump would have
had to reimplement `build_amp2` in Python — which is the very thing the
reference is supposed to check independently.

### 8.2 What moved inside `tagged.json`, and what did not

Rewritten: `channels.li6_alpha.model` and `channels.deuteron.model` — i.e.
`n_of_kc`, `struck_populations`, `p2_moment`, `norm`,
`population_integrated`, `vector_dilution`, `tensor_dilution`,
`p2_moment_mixture_uniform`, plus the unchanged `grid` / `k_pts` / `c_pts`
carried over from the file itself so the SAMPLING never moved, only values.
Worst move against the old file: **+6.277e+02 relative** (⁶Li,
`struck_populations`) and **+8.905e+02** (deuteron) — probabilities near zero
in the old file, so the relative number is large; the absolute
`struck_populations` move is ≤ 0.86.

**Untouched, verified equal to the byte against `HEAD`:** the whole of
`channels.li7_alpha.model` (§3 — one L = 1 wave, nothing moved, so the port
gate against polligen is still good and the re-pin script **re-checks** it at
rtol 1e-12 rather than overwriting it; measured worst 9.591e−13, on
`p2_moment_mixture_uniform`, a quadrature the test compares at 1e−9), and for
**all three** channels `waves`, `base`, `beam_configs`, `boost_spectator`,
`label`, `j_ion`, `s_struck`, `s_spec`, `s_channel`, plus the top-level
`P_D_LI6` and `P_D_DEUTERON`. The only new top-level keys are `provenance`
and `provenance_note`.

The rtol-1e-12 gate (`tests/test_tagged.cpp`, *"the model grid and its tables
against polligen"*) was re-run against the new file and **passes** — including
`norm` at the README's looser 1e-9, which the move would otherwise have broken.
*(That move was published here as "+2.0e−05"; re-measured 2026-09-06 as the
re-pinned file against `git show HEAD:validation/reference/tagged.json` it is up
to **5.81e−05 relative** — ⁶Li M = 0, 1.000026689626 → 0.999968571520 — and
3.97e−05 on the deuteron. `population_integrated` likewise moves up to
**3.02e−06 absolute**, not the ≤ 5e−07 the investigation's §6.1 quoted, and
`p2_moment_mixture_uniform` moves on every spin-1 channel (⁶Li Hulthén
−5.3160743e−05 → −5.2153460e−05) where that table recorded a zero. All four
were re-pinned; the corrected magnitudes are in `validation/README.md` and in
`07_cw_sign_investigation.md` §6.1's correction box.)*

### 8.3 Every other reference JSON is byte-identical

SHA-256, before the fix and after the re-pin:

| file | before | after |
|---|---|---|
| `b1_default_li6.json` | `d7bd8ce6a01b581d…` | **same** |
| `beams.json` | `00d24891685d9fe6…` | **same** |
| `bookkeeping.json` | `2db57d2168fc8fcc…` | **same** |
| `coherent.json` | `f77dcf04681e11f7…` | **same** |
| `spectator.json` | `167f1a7549bf6fce…` | **same** |
| `spin.json` | `60a545b10ee3254d…` | **same** |
| `xsec.json` | `05d6d66c847c1f9b…` | **same** |
| `tagged.json` | `046275b4238ce5b1…` | `63ccf2092d0a06e0…` **(the re-pin)** |
| `_manifest.json` | `a0ff5177a448187e…` | `e27cb830b0774446…` **(one line: the generator of `tagged.json`)** |

---

## 9. The four tests the investigation named, and what each became

| test | what it was | what it is now |
|---|---|---|
| *"the model grid and its tables against polligen"* | rtol-1e-12 port gate against polligen for all three channels | same gate, annotated: for ⁶Li and the deuteron it now reads the **re-pinned** blocks and is a pin on the fixed C++; for ⁷Li it is still the polligen port gate |
| *"the Cosyn-Weiss deuteron tensor gate (CW TABLE II)"* | ran on **Hulthén**, applied `A_T∥ = −2 A_zz^wf`, and matched CW's Eq. (6.14) **minimum** at f₂/f₀ = 0.7053 while calling it Eq. (6.13)'s **maximum** at √2 | **replaced by the investigation §8 gate**, on the **AV18** control CW actually used: (a) Eq. (6.12) as an identity to 1e−12 on both channels, (b) the P₂ factorization re-pinned to −1.99928, (c) CW's own k landmarks 0.298 / 1.035 GeV, (d) the envelope extremes −2 / +1 with `env_min × (−1/2) = +1`, (e) TABLE II row by row at the cell centres, (f) the stated range `[−2, 1]`, (g) a **regression guard** that fails if the phase is ever dropped again |
| *"the deuteron S/D interference shape"* | asserted an **interior peak**, `0.2 < k[argmax\|A_zz\|] < 0.45` | re-pinned: the corrected \|A_zz\| is **monotone** on this channel — measured 279/279 non-decreasing steps, worst drop `0.000e+00` — so the maximum sits at the top cell k = 1.2 at **+0.996607**, saturating at CW's +1 as f₂/f₀ → 1.2866. The threshold suppression (0.006418 < 0.02) and the CW-window O(1) statement survive; the window is now pinned as an interval, [0.731102, 0.957798] |
| *"the acceptance-weighted curve reduces to the cell curve"* | `ref[jk] < −0.4` and `weighted[jk] > +0.4` | **both signs re-pinned**: `ref[jk]` **+0.896640**, `weighted[jk]` **−0.953283** (was −0.491506 / +0.522555). The STATEMENT is unchanged and is the one worth keeping — the YR optics sculpt the sample into the longitudinal half of the sphere, where P₂ has the opposite sign to the θ_k = 90° cell, so the weighted curve is the cell curve's mirror. Only the polarity of the pair moved |

Every re-pinned value carries a comment naming this fix and the pre-fix
number. A fifth site, T25's four dilutions, is recorded in §7.2.

> **ONE TOLERANCE WAS WIDENED, and it was not recorded here until 2026-09-06
> (verification pass).** In the rewritten Cosyn–Weiss gate, TABLE II's
> θ_k ≈ 90° cell-centre row went from `CHECK_CLOSE_AT(a_par_90, 0.9997, 0,
> **1e-3**)` in the old gate to `CHECK_CLOSE_AT(azz_tensor_curve(v, ic90)[ik30],
> +0.99931, 0, **2e-3**)` in the new one — the only loosened tolerance in the
> diff. **The reason, measured:** the pin is now the cell value itself and the
> residual against it is **2.77e−06** (measured 0.9993127695 against the pinned
> +0.99931), so 2e−3 is slack by ×722 and the old 1e−3 would still have passed
> by ×361; what makes a literal tolerance here awkward is CW's own **+1**,
> which the cell centre misses by **6.87e−04** — within a factor 1.5 of the
> 1e−3 the old row used. The tolerance is not load-bearing either way; the
> load-bearing assertion for this row is the identity (a), residual 8.882e−16,
> and the regression guard (g). **Everything else in the diff is equal or
> tighter**: 0.02 → 0.01 on the √2 landmark, 3e−3 → 1e−3 on the factorised
> extremes, T25 stays at rtol 1e−6, and the port gate's rtol 1e−12 / `norm`
> 1e−9 are untouched. Of the other two TABLE II rows, the θ_k ≈ 0 one was
> already at 2e−3 before the rewrite (`a_par_0` against −1.9378) and the
> k ≈ 1.00 GeV one is new.

Tests that keep passing untouched, as the investigation predicted: the pure-S
and pure-D limits, the whole ⁷Li block, *"every channel is normalized"*,
*"sum_M n_M is isotropic"*, `boost_spectator`, the dilution tests, T24 and
T26–T27, *"the VMC S-D interference flips sign below the S node"*,
*"sample_kc reproduces its own density"*, and
`test_pipeline.cpp`'s *"A_zz^tag(k) closes on the acceptance-weighted truth"*
(a generator-vs-own-model closure, sign-blind).

---

## 10. Documentation sites rewritten

| site | before | after |
|---|---|---|
| `docs/USAGE.md` §"⁶Li A_zz^tag at k = 0.20 GeV" row | `+0.845 \| +0.452` | **`−1.207 \| −0.519`**, with a correction box giving the measured cell (k = 0.1979), the optics, and the statement that every other row of that table is unchanged and the spin-blind rate identical to 15 digits |
| `docs/CONVENTIONS.md` §VMC sign paragraph | *"negating the D table flips A_zz^tag from +0.45 to −0.52"* | the sentence survives as an interval (**"flips A_zz^tag by 0.97"**) with the polarity corrected in place: shipped **−0.5191**, negated **+0.4518** — the two numbers swapped roles |
| `docs/PHYSICS_CHANNELS.md` A_zz^tag row | reconciliation values `0.845 / 0.511 / 0.452 / −0.519` | **`−1.207 / −0.600 / −0.519 / +0.452`**, stating that the old `vmc-flipD` control column is now the physics and vice versa. The row's *"VMC and Hulthén agree in sign inside the accepted window"* reading **survives** — both are negative there |
| `docs/open_items/vmc_reconciliation.md` §"Tagged tensor asymmetry" | the 2026-08-29 table | **the table STAYS as measured**, with a dated **correction box** below it carrying the re-measured five rows and naming the `vmc-flipD` ↔ physics swap. A dated record is not overwritten |
| `validation/README.md` §`tagged.json` | *"dumped from polligen"* | the full re-pin record of §8, generated from `dump_polligen_reference.py`'s own `doc()` text so it survives a regeneration |

**Second pass, 2026-09-06 — every remaining site that quotes an inheriting
number, live or dated.** A live site gets the corrected number with the fix
named and dated; a dated record keeps its original and gets a correction box
beside it. Nothing was silently overwritten.

| site | kind | what changed |
|---|---|---|
| `docs/OPEN_ITEMS_SOLUTIONS.md` §1 (two bullets) | live | *"⁶Li A_zz^tag at k = 0.20 GeV +0.845 → +0.452"* → **−1.207 → −0.519**, with the cell, optics and `n_phi` named. And *"Tagged A_zz^tag(k) roughly halves"* is **withdrawn as a reading**: re-measured, the VMC/Hulthén ratio is 0.430 / 0.273 / 0.209 / 0.113 / −0.200 over the five k (it was 0.535 / 0.398 / 0.348 / 0.234 / −0.554), i.e. a factor **2.3 at k = 0.20 falling to 8.8 at k = 0.40** — not a halving |
| `docs/OPEN_ITEMS_SOLUTIONS.md` §11.4 | live | a box measuring **all five C5 bullets** against the fix and finding none of them inherits it — each quantity there is angle-integrated or a `constexpr` closed form in P_D. Numbers in the box, not adjectives |
| `README.md` "External anchors reproduced" | live | *"Cosyn–Weiss deuteron TABLE II (+1/−2)"* → the anchor named properly, with the **mapping corrected to +1** and the identity residual **8.881784e−16 over 26 880 cells** against 2.740499e+00. Also records the one deliberate exception to the *"agree with polligen at rtol 1e-12"* line above it |
| `docs/DEVELOPMENT_PLAN.md` §4 row 4 | live | the tagged gate's description rewritten to what it now checks (AV18 control, Eq. 6.12 as an identity, CW's own k landmarks 0.298121 / 1.034872 GeV), with the old *"peak k = 0.31 GeV, A_T∥ extremes +1/−2"* named as the Hulthén-channel error it was |
| `python/lipolgen/cli.py` `--cluster-wave` help | live, **user-facing** | a SIGN NOTE: any `A_zz^tag` from a run before 2026-09-06 has the wrong sign; the two current values; and the explicit list of what is *not* affected (rates, tag fractions, P_D, both dilutions, all of ⁷Li) |
| `python/lipolgen/configs.py` module docstring | live | records that this module **does not** inherit the fix and why — its `tensor_dilution` is `ClusterConfigSampler`'s closed form `1 − 0.9 P_D`, a different quantity from the tagged model's grid integral, and phase-blind |
| `docs/benchmarking/BENCHMARK_PLAN.md` row 2 and §5 item 1 | dated plan | marked **RESOLVED AND APPLIED**, with every clause of the item confirmed and the measured outcome |
| `docs/benchmarking/00_in_tree_checks.md` gate **E-1** | dated survey | correction box: criterion (c) was wrong twice over; the three pass values (0.99940, k = 0.3098, 0.9997/−2) are pre-fix; what the rewritten gate measures instead |
| `docs/benchmarking/06_critic.md` **U-2**, `04_theory.md` **T-6** | dated survey | the refutation is annotated **repaired**; T-6 gains the settled mapping and the note that the gate moved onto the AV18 control |
| `docs/surveys/needs_survey.md` §5.2, §5's item 5, and the anchor table | dated survey **of the sibling** | three correction boxes. The §5.2 entry told a replacement generator to reproduce `test_cosyn_weiss_tensor_gate` — the gate that is wrong twice over — so it now says **do not**, and states what to reproduce instead |
| `docs/open_items/run_2026-09-03/phase_C_numbers.md` §C5.2, §C5.4 | dated record | §C5.4's *"the relative S–D sign does not flip"* is the **third** site carrying the false step (`07` §5 named two). Box separates what stays true (ψ₂ = +W, every clause) from the unstated corollary that was false. §C5.4's dilution table and §C5.2's MC-error column both re-measured and **unmoved** |
| `docs/open_items/run_2026-09-03/AUTHOR_DECISIONS.md` rows 9/10/11, §B9, §B10 | dated registry | **row 9's evidence moved** — its (c) leg *"the control's job is to be the Cosyn–Weiss tagged limit"* now argues the other way, since the CW gate runs on AV18; and *"`tagged.json` pins Hulthén"* is weaker now that the spin-1 blocks are self-generated. **No decision taken.** Rows 10 and 11 checked and annotated: neither inherits |
| `docs/open_items/run_2026-09-03/SUMMARY.md` | dated record | the +2.069 % / −2.027 % bullet re-measured (unmoved), and a new finding **6b** recording that a shipped tensor number was wrong and is now fixed |
| `docs/USAGE.md`, `docs/open_items/vmc_reconciliation.md`, this file §7.1 | live + dated | the spin-blind-rate claim tightened from *"identical to 15 digits"* to the bound it actually supports — **≤ 4.0e−16 relative, 1–2 ulp, VMC bit-identical** — after an independent re-measurement found the VMC fraction's 15th digit is summation-order dependent (`0.033810227625845` as first published, `0.033810227625842` re-measured; the latter reproduces the unnormalized `6.005768857484e+01` exactly) |
| `docs/open_items/run_2026-09-06/phase_CW_numbers.md` §11 | this file | rewritten as the full sibling-impact note: eight inheriting sites and seven checked-safe ones, each with file:line and the measured number |

**Checked and found to carry no inheriting number** (recorded so the search is
not repeated): `docs/open_items/run_2026-09-02/*.md` — `PLAN.md:60, :104` name
the observable without a value; `design_C_tensor_rc.md:877, :2068` are
structural definitions and a self-consistency identity (**T2**, sign-blind:
it compares the generator against its own `azz_tensor_curve`);
`design_D_b1_li6.md:293, :342, :375` are the CG closed form
`1 − (9/10) P_D` and the b₁ᵈ sector, which is `b1_nuclear` and was already
right. The CLI's own `--help` output carries no A_zz^tag number anywhere else
(checked by running it). `python/lipolgen/configs.py` carries none.

---

## 11. What in `PolarizedLithiumSim` inherits the bug

**`PolarizedLithiumSim` was NOT modified — nothing in this section was written
to that tree.** Every line below is a read, and every LiPolGen number beside it
was measured on 2026-09-06 with the before/after instrument of §0: on the fixed
library, negating the stored L = 2 table reproduces the pre-fix amplitude
exactly (it cancels the new `(-1)^floor(L/2)`), which is why every "before"
column here reproduces this tree's own pre-fix published values cell for cell.

### 11.1 The root: the same missing `i^L`, one line

`evgen/polligen/tagged.py:247`, inside `TaggedModel._amp2_table`:

```python
                amp = amp + (self._rad[w.l][:, None] * cg
                             * theta_lm(w.l, int(round(m_l)), cc))
```

No phase — exactly what `src/core/tagged.cpp` `build_amp2` did until
2026-09-06. The fix there is the same one line: multiply by
`(-1)**((w.l // 2) % 2)`. **Parity fixes L mod 2 inside one channel**, so this
is the whole of it; `+1` on L = 0, `+1` on L = 1 (⁷Li: one wave, the common
`i` is a global phase), `−1` on L = 2.

**Scope: the two spin-1 channels only** — `li6_alpha_channel` and
`deuteron_channel`. `li7_alpha_channel` is a single L = 1 wave with nothing to
interfere with and is **exactly** unaffected (measured on this side: ⁷Li's
`p2_moment` is −0.1999349117 / +0.1998480944 before and after, and
`validation/reference/tagged.json`'s ⁷Li `model` block came through the re-pin
**byte-identical**, §8.2).

### 11.2 What carries the inverted sign, with the number

| # | site | what it carries | measured here, before → after |
|---|---|---|---|
| 1 | `evgen/polligen/tagged.py:247` | the missing phase itself | — |
| 2 | `evgen/polligen/tagged.py` `n_of_kc`, `struck_populations`, `p2_moment`, `p2_moment_mixture`, `sample_kc`, `azz_tensor_curve`, `azz_tensor_curve_weighted` | every **angle-differential** tagged quantity, at full size | ⁶Li `p2_moment(M=±1)` **+0.0448133305 → −0.0622510957** (×1.389); `(M=0)` **−0.0897861433 → +0.1243457311** |
| 3 | `evgen/scripts/money_tagged_azz.py:151, 238, 242, 246, 278, 290, 318` (`azz_wf_curve = tagged.azz_tensor_curve`) | **the published money plot**, both panels: the analytic θ_k = 90° curve, the acceptance-weighted truth overlay, and the accepted-sample points, which come from `sample_kc`'s own (M, m_S, cos θ_k) correlation | ⁶Li Hulthén 90° curve at k = 0.3313 GeV: **−0.4786 → +0.9279** |
| 4 | `evgen/money_tagged_azz_6Li.png` | the published figure produced by #3 | inherits #3 in full |
| 5 | `evgen/README.md:113-115` | *"the folded A_zz reads **+0.49** and **−0.07** at k ≈ 0.33 GeV/c where the 90° curve says **−0.48**"* | the −0.48 reproduces here as **−0.4786** pre-fix and is **+0.9279** post-fix. All three numbers in that sentence flip; the *reading* they support — that the accepted sample is not at θ_k = 90° and the swing is θ_k sculpting, not the wave function — is **sign-blind and survives** |
| 6 | `evgen/tests/test_tagged.py:737-790` `test_cosyn_weiss_tensor_gate` | **wrong twice over**, and the reason the bug survived review on both sides: it applies `A_T∥ = −2 A_zz^wf` (docstring `:751-753`, asserts `:784-788`) when the mapping is **+1**, and it runs on `deuteron_channel()` — the **Hulthén** pair, whose f₂/f₀ never reaches √2 anywhere on the grid — while quoting CW's AV18 landmark | its four pinned values are all pre-fix: `ratios.mean() == 0.99940` (`:768`), `envelope[j] == 1.0` at `m.k[j] == 0.3098` (`:775-777`), `a_par_90 == 0.9997` and `a_par_0 == −1.9378` (`:785-786`). LiPolGen's replacement checks Eq. (6.12) as an **identity** on AV18: **8.881784e−16** over 26 880 cells, against **2.740499e+00** |
| 7 | `evgen/tests/test_tagged.py:233-244` `test_deuteron_tagged_azz_wf_shape` | pins the pre-fix **shape**; would fail against a fixed implementation | `assert 0.2 < peak < 0.45` (`:243`) — the peak of \|A_zz^wf\| at θ ≈ 90° moves **k = 0.3098 → 1.2000 GeV** (post-fix \|A_zz\| is monotone in k on this channel: **279/279** non-decreasing steps, worst drop 0.000e+00), so it **breaks**. `assert …max() > 0.45` (`:242`) survives, but on **0.499835 → 0.957798**. The threshold clause `abs(a[k<0.02].max()) < 0.02` (`:240`) survives at **0.006397 → 0.006418** |
| 8 | `evgen/tests/test_tagged.py:93-136` `test_density_matches_independent_transcription` | an independent re-transcription of the same formula — so it re-transcribes the same omission and agrees with it | would need the phase added on both sides |

### 11.3 What is in that tree and does **not** inherit — checked, not assumed

| site | why it is safe | measured |
|---|---|---|
| `fastsim/polli_fastsim/spectator.py` (whole module) | **spin-blind by construction.** It samples \|ψ_L(k)\|² one partial wave at a time (`momentum_density(k, kappa, beta, l_wave)`, `:212`; `sample_k`, `:222`) and boosts; there is **no ion projection M, no Clebsch–Gordan and no interference term anywhere in the file** — its public surface is `spectator_lab_kinematics` (`:280`) and `breakup_lab_kinematics` (`:317`), both pure kinematics. A single-L density has no relative phase to get wrong | — |
| `fastsim/polli_fastsim/polarized.py:873` `LI6_B1_RANK2_TRANSFER = 0.921947` | quotes the tagged ⁶Li `tensor_dilution`, which is **angle-integrated**: `∫Θ₀Θ₂ dc = 0` by L-orthogonality kills the flipped cross term and only the 96-cell midpoint-quadrature residual survives | **0.9219467089 → 0.9219489770**, **+2.46e−6 relative**, against the **1e−4** tolerance it is pinned at (`evgen/tests/test_tagged.py:725`). **Does not inherit a meaningful error.** The comment at `polarized.py:869-870` does quote the pre-fix digits — *"0.921970 against the 0.9219467 measured"* — which now read **0.9219490** |
| `fastsim/polli_fastsim/polarized.py:883` `b1_li6_from_deuteron`, `LI6_B1_PER_NUCLEON = 2/6`, `LI6_B1_LEGACY_TRANSFER = 0.87` | linear in the transfer above, or closed forms in P_D. P_D is a **norm**, and norms are phase-blind | vector dilution **0.8699393995 → 0.8699431789** (+4.3e−6 relative, pinned at 2e−3) |
| `evgen/scripts/tagged_polarimetry_7li.py`, `evgen/tagged_polarimetry_7Li.png` | ⁷Li only — one L = 1 wave | **exactly** unchanged |
| `evgen/scripts/target_mass_bound.py`, `fastsim/polli_fastsim/farforward.py`, `beams.py`, `asymmetries.py`, `structure.py`, `delta_models.py` | no tagged spin amplitude; geometry, constants, inclusive structure functions | — |
| `plans/01_findings_physics_case.md:249` | *"Tagged tensor asymmetries O(1) at ~300 MeV/c"* — a claim about **Cosyn–Weiss's paper**, not about polligen's output, and it states a magnitude with no sign | survives; the fix moves ⁶Li A_zz^tag(0.20 GeV) from 0.845 to **1.207** in magnitude, so "O(1)" is if anything better supported |
| `plans/02_phase1_event_generation.md:32, 36, 187, 218, 251, 275`; `plans/07_plb_letter_gluonometry.md:21, 41, 220, 251` | every `A_zz` / `δA_zz` there is the **inclusive b₁** asymmetry (`−(2/3) b₁/F₁`), a different observable computed on a different path; `plans/07:41` explicitly **excludes** tagged A_zz from that letter | — |

### 11.4 Reports — qualitative, and they survive

The published reports carry the tagged tensor asymmetry as a *magnitude* claim
and a figure reference, never as a signed number, so none of them states
something now false:

- `reports/nanowire_far_forward.template.html:142` (and the built `.html` /
  `.pdf`) — *"the α-tagged tensor asymmetry A_zz^tag(k) … the S/D interference
  of the α–d relative motion"*, routed *"through … `money_tagged_azz.py`"*.
  The prose is right; **the figure it points at (#4 above) is not.**
- `:252` — *"the tagged tensor asymmetry is a one-point measurement at the
  O(1) end of the S/D interference"*. Sign-blind, survives.
- `reports/polarized_li_primer.template.html:593, 633` — *"rides the α–d D
  wave"*, *"the order-unity tagged tensor asymmetries (§5.4)"*. Sign-blind,
  survive.
- `reports/polarized_li_primer.template.html:751` — *"The Cosyn–Weiss
  closed-form tagged tensor asymmetry (§5.4) is added as a **validated** limit
  of the generator."* **This one is now false as written**: the validation it
  refers to is `test_cosyn_weiss_tensor_gate` (#6), which is the gate that was
  wrong twice over. The same sentence's *"⁶Li b₁ transfer is the rank-2 tensor
  dilution ⅓ × 0.92"* is fine (§11.3, row 2).

### 11.5 The consequence for the port gate

Because `polligen`'s `_amp2_table` carries the bug, **`tagged.json` can no
longer be regenerated from it** for the two spin-1 `model` blocks. Those blocks
are now dumped from this library by `validation/repin_tagged_from_lipolgen.py`
and carry `provenance: "LiPolGen post-fix, formerly polligen"`; ⁷Li's block and
everything outside `model` are still polligen's, byte for byte (§8). Both
scripts under `validation/` say so, and so does `validation/README.md`.

**The right fix on that side is the same one line at
`evgen/polligen/tagged.py:247`.** Until it lands, the two trees disagree on
every angle-differential spin-1 tagged quantity, and **this tree is the one
that agrees with Cosyn–Weiss Eq. (6.12)** — as an identity, to 8.88e−16.

---

## 12. Verification pass, 2026-09-06 — what three adversarial lenses left behind

Three lenses (*only-what-should-move*, *gate-quality*, *propagation*) were run
against the fix and the documentation pass above. **None refuted the code fix**:
the phase is the specified one, only the two spin-1 tagged channels move, only
`tagged.json` moved, and Cosyn–Weiss Eq. (6.12) holds as an identity. What they
refuted were CLAIMS around it. Every item below is fixed at its site, and every
number in this section was re-measured in this pass.

| what was wrong | where it was | what it now says |
|---|---|---|
| *"tag fractions are unaffected"* | `python/lipolgen/cli.py` `--cluster-wave` help (user-facing), §7.1 and §7.2 above, `07` §6.1/§6.5, `run_2026-09-06/STATUS.md` row 9 | A tag fraction is an acceptance-weighted integral over n_M: invariant only for a spin-blind or category-averaged fill. Uniform-M mix unmoved to **15 digits** (0.024675932148828 / 0.033810227625842), equal-thirds average +3.1e−9; but `tensor-thirds` categories **0.028127 → 0.018403** (azz±) and **0.017775 → 0.037222** (azz0), and the CLI's own `--pz 0.7` ladder **0.027030 → 0.020396** (−24.5 %). ⁶Li VMC 0.035296 → 0.032155, d Hulthén 0.025802 → 0.021492, d AV18 0.022175 → 0.019407; ⁷Li exactly unchanged |
| *"`USAGE.md`'s tag-fraction rows … were left alone"* (§7.1) | `docs/USAGE.md` (two tables), `docs/open_items/vmc_reconciliation.md` | Both **REGENERATED on the fixed build** by their own documented command (`vmc_tag_fractions.py --events 40000 --configs 0,1,2`, seed 20260829): ⁶Li YR high-acceptance **0.0301 / 0.0264 / 0.0279**, ⁶Li tagging **0.3451 / 0.2551 / 0.3145** (were 0.0286/0.0249/0.0266 and 0.3410/0.2530/0.3115); every ⁷Li row identical. The captions say the build and carry the invariance statement |
| the false step, still live in two primary docs | `docs/USAGE.md`, `docs/PHYSICS_CHANNELS.md`, `docs/OPEN_ITEMS_SOLUTIONS.md` §C5.4 bullet, `run_2026-09-03/STATUS.md` row 9, plus the two test comments that assert the stored convention | ψ₂ = +W is the **stored** convention and stays; `build_amp2` consumes φ₂ = i²ψ₂ = **−W** and applied no phase until 2026-09-06 |
| stale digits inheriting the quadrature residual | `USAGE.md` deuteron-control dilutions; `CONVENTIONS.md:93`, `PHYSICS_CHANNELS.md`, `beams.hpp`, `tagged.hpp`, `test_tagged.cpp` (comment) for the ⁶Li vector dilution; `constants.hpp` and `design_G_cluster_config.md` for 0.9219467; `breakup.hpp` for 0.932495 / 0.913595 | 0.932496109 / 0.959488878 / 0.913595979 / 0.948146339; **0.8699431789** against the closed form 0.869950, **7.84e−6** (was 1.22e−5); **0.9219490** (was 0.9219467); 0.932496 / 0.913596 |
| the reconciliation record was one script run from losing its correction | `validation/vmc_reconcile.py` (template sentence, polarity inverted) and `validation/vmc_tag_fractions.py` (rewrites everything after the APPEND marker, where the hand-inserted box lived) | Both **GENERATORS** carry it now: the sentence reads −0.52 → +0.45 with the measured cells, the correction record and the tag-fraction invariance note are emitted by `vmc_tag_fractions.py`, the `vmc-flipD` column is relabelled *(= the PRE-FIX physics column)*, and the box's `:131-133` citation is replaced by a section anchor plus the correct `:133-136`. Regenerated and verified **idempotent** |
| the re-pin record's magnitudes | `validation/README.md` and its generator `dump_polligen_reference.py`; `07` §6.1; §8.2 above | `norm` up to **5.81e−05 rel** (was "+2e-5"), `population_integrated` up to **3.02e−06 abs** (was "≤5e-7"), and `p2_moment_mixture_uniform` — **absent from the list**, recorded as a zero in `07` §6.1 — moves on every spin-1 channel (⁶Li Hulthén −5.3160743e−05 → −5.2153460e−05, 1.9e−2 rel) and was re-pinned. Measured as the re-pinned file against `git show HEAD:validation/reference/tagged.json` |
| *"no tolerance was widened"* | §9 above | One was: the CW gate's TABLE II θ_k ≈ 90° row, atol **1e−3 → 2e−3**. Recorded with its measured residual (2.77e−6, slack by ×722) in §9's box; the test was not touched |
| registry board not updated | `run_2026-09-03/STATUS.md` rows 9 and 10 | Both cells amended to §B9(c)'s "evidence moved" box and §B6(c)'s re-measured residual; no decision taken, both rows stay open |

**Three more live sites carried the pre-fix tag-fraction SAMPLE** and were
brought onto the regenerated numbers with the reason named: `README.md`'s
`--cluster-wave` row (⁶Li α-tag fraction 0.0249 → **0.0264** at 10 × 99.5),
`PHYSICS_CHANNELS.md`'s VMC-tables row (0.0348 vs **0.0264**, ×**1.315**, and
0.2486 vs **0.2551** at the tagging optics — it said ×1.40 against 0.0249) and
`OPEN_ITEMS_SOLUTIONS.md` §1 (three bullets, including the two-configuration
check 0.0301 → 0.0365 at 5×41 and 0.0279 → 0.0348 at 18×275). One unrelated
anchor slip found while checking them and fixed: `PHYSICS_CHANNELS.md`'s ⟨P₂⟩ =
−T/5 row cited `README.md:109` (anchor read 2026-09-15 "rtol 1e-12"), which is the rtol-1e-12 sentence; the statement
is at `README.md:115` (anchor read 2026-09-15 "2025 Table 1").

**Not changed by this pass:** `src/`, `validation/reference/` (SHA-256
unchanged; `repin_tagged_from_lipolgen.py --check` still reports worst move
0.000e+00 for both re-pinned blocks and 9.591e−13 for ⁷Li), the CW gate itself,
and `PolarizedLithiumSim`. Edits under `include/` and `tests/` are
**comments only** (verified by a non-comment diff).

**End state, re-measured after the pass:** C++ **404 cases / 404 passed / 0
failed / 1 skipped, 17 240 391 assertions**; pytest **1004 passed / 127
skipped**; docs gate **1220 references (strict) / 96 ranges / 7 external / 0
broken / 6 allow-listed** (+19/115 `SPIN32_FINITE_GAMMA.md`, +8 external
`PYTHIA_BRIDGE.md`); SPDX **102/102**.

**How the gate was brought back to 0 broken, and what that cost.** The
correction boxes shifted 127 citations. `--fix` re-anchored **121** of them
(pure line shifts). The rest were re-pointed **by hand and by CONTENT**, each
verified against its recorded `sha256` before the citation was moved:
`USAGE.md`'s four cited blocks (§4 opt-in table, §5 far-forward, §7
checks-for-free, §8 command-line generators — `--fix` will not re-anchor a
range whose recorded start line has become blank), `breakup.hpp:161 → :164`
(the ³He-control comment), `cli.py:1354 → :1368` (`rc_model`, a pinned use
site) and `cli.py:1117 → :1131`. One citation was *extended* rather than
moved: `include/lipolgen/tagged.hpp:54-60 → :54-62`, because the comment it
cites is two lines longer now and the old range ended mid-sentence.

Against the sidecar as it stood before this pass, `validation/physics_channels_ranges.json`
holds the same **219** fingerprints, of which **14 moved with byte-identical
content** and exactly **two carry new content**, both deliberate and reviewed:
`include/lipolgen/tagged.hpp:54-62` (the re-measured residual in the "one wave
function" rule) and `docs/open_items/vmc_reconciliation.md:197-206` (the
regenerated A_zz table, cited from `PHYSICS_CHANNELS.md`'s A_zz^tag row, which
now says what that table is). **No cited block was edited to make the gate
pass.**
