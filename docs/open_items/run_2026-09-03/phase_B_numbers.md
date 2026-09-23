# Phase B — every number, measured

Task B of `docs/open_items/run_2026-09-03/PLAN.md`. Written after the runs, from
output quoted with the command that produced it. Nothing here is an estimate;
where a number is *not* measured the sentence says so.

Baseline at the start of Phase B: commit `ac22331`, tree clean, **374 doctest
cases / 17 208 032 assertions / 1 skipped / 0 failed**, **261 pytest passed**,
**1049 doc references / 0 broken**.

---

# B1 — the ⁶Li elastic C0 shape is a BAND, and the sign it decides

## B1.0 What this task did, and what it deliberately did not

`include/lipolgen/rc.hpp`'s `LI6_FF_HO_A_FM = 1.9069` / `LI6_FF_HO_ALPHA =
0.13822` are, by their own comment, "starting parameters … TO BE REFIT", and
`run_2026-09-02/phase_C_numbers.md` §8.3 records the refit as **NOT DONE, Q1
stays open**. That is still true after this task. **Q1 is not closed here.**

What is new is a **price tag**: a second, independent C0 shape, run as the other
edge of a band, and the measurement of what the shipped default's unfitted shape
was silently deciding. It was deciding a published sign.

**Why a band and not the refit.** A two-parameter harmonic oscillator cannot be
refitted to the second shape while holding ⟨r²⟩ — its high-q fall is
exp(−q²a²/4) with `a` fixed by ⟨r²⟩, so the only way to lift q ≈ 2.5 fm⁻¹ by the
required factor 6 is to shrink `a` and break q → 0. A refit therefore means
inventing a 3–4 parameter functional form, a new provenance argument for it, and
a rewrite of T11's closed-form assertions — against elastic data that are still
not in this repository in any machine-readable form. Running two edges prices the
same uncertainty and asserts nothing new.

## B1.1 The new knob

`C0Shape { Ho, VmcFt }` on `HoSpin1FFOptions` (`include/lipolgen/rc.hpp:475`,
field at `:321`) and on `RcOptions` (`:790`), read at exactly one place,
`c0_point` (`src/core/rc.cpp:1061`), which both `HoSpin1FF::fc` and
`HoSpin1FF::fq` call — they share **one** monopole, so the knob moves both
together and cannot be rescaled out of one run. Plumbed through the same
five-file chain `fq_scale` occupies: `RcOptions` → `src/core/rc.cpp` (`RcModel`'s
default form factor) → `python/bindings.cpp` (the `C0Shape` enum, `c0_shape_name`,
both `def_readwrite`s and `meta["rc_c0_shape"]`) → `python/lipolgen/cli.py`
(`--rc-c0-shape {ho,vmc-ft}`) → `python/lipolgen/__init__.py`
(`RC_C0_SHAPES`, `make_config(rc_c0_shape=…)`). Gated by
`python/tests/test_rc.py::test_c0_shape_band_shares_its_anchors_and_flips_the_tail_sign`
and `…::test_c0_shape_reaches_the_tail_and_the_npz_meta`, and by T11's new
`C0Shape::VmcFt` subcase in `tests/test_rc.cpp`.

**`Ho` is the default and nothing about the shipped run moved.** Measured, on a
20 000-event inclusive run at seed 20260713, `--rc tensor-band` with the two
edges: `x`, `Q²`, `y`, `Event::weight`, `rc_tensor_lo` and `rc_tensor_hi` are
**bit-identical**; only `rc_tail` differs (mean 1.020848536828 → 1.020934589099).
The `--rc off` npz is untouched — the whole rc `meta` block, `rc_c0_shape`
included, is written only when RC is on.

## B1.2 What `VmcFt` is

F_point(q) = Σ ρ_p(r) j₀(q·s·r) r² / Σ ρ_p(r) r² over
`data/vmc/density/li6.density` (ANL AV18+UX point-**proton** density, 200 k VMC
samples, R = 0.05…20.05 fm step 0.1, so the plain sum **is** the midpoint rule on
[0, 20.1] fm and the equal bin width cancels in the ratio). Normalised to
F(0) = 1 by construction. The radial rescale

    s = sqrt(LI6_R2_POINT_FM2 / ⟨r²⟩_table) = 1.0089589

is applied so the shape carries the **measured** ⟨r²⟩_point = 6.0788 fm² exactly.

| quantity | measured |
|---|---|
| rows read (`read_anl_plain`), of which non-zero | 201, **99** |
| `4π∫ρ r² dr` recomputed from the file | **2.99922** against Z = 3 (0.03 % low) |
| ⟨r²⟩ of the file itself | **5.971328 fm²**, r = **2.44363 fm** |
| the file's own printed rms (`LI6_R_POINT_VMC_FM`) | 2.4433 fm |
| design target r_point = √6.0788 | 2.4655 fm — the file misses it by **0.9 %** |
| rescale actually applied | **1.0089589** (from the quadrature's own ⟨r²⟩, *not* from the printed 2.4433, which would give 1.0090952 and land ⟨r²⟩ 2.7e−4 off — T11 gates it at 1e−4) |

**The anchors are identical on both edges**, which is the point: the band is on
the unfitted *shape* alone and never on a measured normalisation.

| anchor | `ho` | `vmc-ft` |
|---|---|---|
| F_c(0) | 3 | 3 |
| F_m(0) | 4.907652464 | 4.907652464 |
| F_q(0) | −65.91414947 | −65.91414947 |
| ⟨r²⟩_point, `−6 dF/dq²\|₀` on the unfolded shape | 6.078840 fm² | 6.078796 fm² |

## B1.3 How far apart the two edges are

`F_point(q)`, unfolded, F(0) = 1:

| q [fm⁻¹] | `ho` | `vmc-ft` | vmc/ho |
|---|---|---|---|
| 0.5 | 7.7598e−01 | 7.8075e−01 | 1.006 |
| 1.0 | 3.6097e−01 | 3.9385e−01 | 1.091 |
| 1.5 | 9.9044e−02 | 1.4285e−01 | 1.442 |
| **2.0** | 1.5381e−02 | 3.8576e−02 | **2.508** (×6.29 in F²) |
| **2.5** | 1.1912e−03 | 7.3263e−03 | **6.150** (×37.8 in F²) |
| 3.0 | 1.7719e−05 | 1.2812e−03 | 72.3 |
| 3.0999 (`ho`'s C0 zero) | −1.3368e−08 | 9.1178e−04 | — |
| 3.5 | −4.0072e−06 | 2.3348e−04 | — |
| 10 | −3.1131e−39 | 5.7106e−14 | — |
| 88 | −0 (underflow) | 2.6167e−129 | — |

They agree to 0.6 % at q = 0.5 fm⁻¹, where both are pinned to the same ⟨r²⟩, and
diverge exactly where the elastic tail lives (q = 2 fm⁻¹ ⇔ t ≈ 0.155 GeV²).
**The `ho` edge changes sign at q₀ = 3.0999 fm⁻¹ and the `vmc-ft` edge never
does.** That difference, not the size of the shape, is what moves σ^el_T.

## B1.4 The high-q continuation — its own model choice

The t-peak integrand runs to η_A = S_x/4M_A², i.e.

| x | y at Q² = 5 | η_A range | q range [fm⁻¹] |
|---|---|---|---|
| 0.01 | 0.12563 | 6.96e−07 … 2.39e+01 | 0.047 … **277.6** |
| 0.03 | 0.04188 | 6.28e−06 … 7.97e+00 | 0.142 … 160.2 |
| 0.10 | 0.01256 | 7.04e−05 … 2.39e+00 | 0.476 … **87.8** |
| 0.30 | 0.00419 | 6.38e−04 … 7.97e−01 | 1.43 … 50.7 |

Nobody has data there and the VMC table is noise long before it. Propagating the
file's own `DRHORP` column (bins treated as independent) through the transform:

| q [fm⁻¹] | F | ±1σ (MC) | F/σ |
|---|---|---|---|
| 2.0 | 4.06394e−02 | 2.46e−04 | **165** |
| 2.5 | 7.94685e−03 | 1.99e−04 | 40.0 |
| 2.82 | 2.45446e−03 | 1.76e−04 | 13.9 |
| **3.0** | 1.37768e−03 | 1.65e−04 | **8.3** |
| 3.0999 | 1.07570e−03 | 1.60e−04 | 6.7 |
| 3.5 | 6.61274e−04 | 1.42e−04 | 4.7 |
| 4.0 | 3.60610e−04 | 1.24e−04 | 2.9 |
| 4.3 | 1.69582e−04 | 1.15e−04 | 1.5 |
| 4.5 | 7.98659e−05 | 1.10e−04 | 0.7 |
| 5.0 | −2.89466e−06 | 9.92e−05 | −0.03 |

(unrescaled q; the text-rounding of the density column — 4–5 significant figures
— contributes ≤ 7.4e−07 at every q above, two to three orders below the MC
error, so it is not what limits this.) Above q ≈ 4.3 fm⁻¹ the transform **is**
its own error bar, and the 0.1 fm grid aliases anyway above the Nyquist
q = π/Δr = 31 fm⁻¹.

**The choice made.** Direct transform up to `LI6_VMC_C0_Q_CUT_FM` = 3.0 fm⁻¹ —
the last q at which the table is still signal-dominated — and above it

    F(q) = F(q_cut) · exp(−λ (q − q_cut)),  λ = 3.404842 fm

with λ the transform's own **mean** log-slope over
[`LI6_VMC_C0_Q_MATCH_FM`, `LI6_VMC_C0_Q_CUT_FM`] = [2, 3] fm⁻¹, the last decade
over which it is signal-dominated. Continuous in **value**, deliberately **not**
in slope: the *local* log-slope at 3.0 fm⁻¹ is already 2.645 fm and at 3.5 fm⁻¹
it is 0.77 fm, which is the noise flattening the tail and not the physics.

Two consequences, stated rather than discovered downstream: the continuation is
positive and monotone, so **the `vmc-ft` edge has no C0 zero at any q**; and the
continued tail is numerically dead (F = 5.7e−14 at q = 10 fm⁻¹), so the
continuation **law** is not what moves the numbers. Measured — the whole table
rerun with a lower cut and with a Gaussian law matched to the same mean
log-slope, against an independent reimplementation of the transform:

| x (Q² = 5) | exp, q_cut = 3.0 (**shipped**) | exp, q_cut = 2.5 | Gaussian, q_cut = 3.0 |
|---|---|---|---|
| 0.01 | −5.470470e−04 | −5.470505e−04 | −5.470470e−04 |
| 0.03 | −8.042604e−04 | −8.042675e−04 | −8.042604e−04 |
| 0.10 | −4.070157e−05 | −4.073873e−05 | −4.070139e−05 |
| 0.30 | +1.555578e−02 | +1.555129e−02 | +1.555581e−02 |

Worst case **0.09 %** (x = 0.10, q_cut 3.0 → 2.5). The continuation is a stated
model choice with a measured, negligible cost — and the independent
reimplementation reproduces the shipped edge to every printed digit, so the two
implementations of the transform agree.

## B1.5 What moves — the republished §8.1c rows

Per-nucleon `RcModel::tail_sigma_at(x, Q²)` and `tail_ratio_at(x, Q², q_N)`
through the shipped Python API (config 1, ⁶Li 10 × 99.5 GeV/u ⇒ s = 3980 GeV²,
`n_eta = 128`, `--plan tensor-thirds --pzz 0.6`). **Every `ho` row reproduces
`run_2026-09-02/phase_C_numbers.md` §8.1c exactly**, which is what licenses the
`vmc-ft` rows beside them.

| x | Q² | edge | `(1/6)σ^el_T/σ^el_U` | `σ^q_U/σ^el_U` (S on) | `σ^el_U` [GeV⁻²/nucleon] | `ΔA_zz` (tail) |
|---|---|---|---|---|---|---|
| 0.01 | 2 | ho | **−5.094305e−04** | 0.28237 | 2.204782e−05 | +8.381444e−08 |
| 0.01 | 2 | vmc-ft | **−5.470471e−04** | 0.27942 | 2.228045e−05 | +7.902630e−08 |
| 0.01 | 5 | ho | **−5.094304e−04** | 0.28237 | 2.221727e−05 | +3.516828e−07 |
| 0.01 | 5 | vmc-ft | **−5.470470e−04** | 0.27942 | 2.245168e−05 | +3.229114e−07 |
| 0.01 | 10 | ho | −5.094304e−04 | 0.28237 | 2.294679e−05 | +1.093833e−06 |
| 0.01 | 10 | vmc-ft | −5.470470e−04 | 0.27942 | 2.318890e−05 | +9.719544e−07 |
| 0.03 | 5 | ho | **−7.358312e−04** | 0.54628 | 1.123928e−06 | — |
| 0.03 | 5 | vmc-ft | **−8.042604e−04** | 0.53452 | 1.148665e−06 | — |
| 0.10 | 2 | ho | **+1.550648e−04** | 2.74787 | 1.331544e−08 | +1.145119e−09 |
| 0.10 | 2 | vmc-ft | **−4.152257e−05** | 2.48840 | 1.470388e−08 | +1.153117e−09 |
| 0.10 | 5 | ho | **+1.559536e−04** | 2.74748 | 1.331742e−08 | +6.675139e−09 |
| 0.10 | 5 | vmc-ft | **−4.070157e−05** | 2.48806 | 1.470596e−08 | +6.716315e−09 |
| 0.10 | 10 | ho | **+1.562196e−04** | 2.74739 | 1.332082e−08 | +2.543549e−08 |
| 0.10 | 10 | vmc-ft | **−4.045848e−05** | 2.48798 | 1.470970e−08 | +2.557536e−08 |
| 0.30 | 5 | ho | +1.328679e−02 | 1000.19931 | 1.336748e−12 | +1.266460e−11 |
| 0.30 | 5 | vmc-ft | +1.555578e−02 | 329.84878 | 4.053416e−12 | +1.733428e−11 |
| 0.30 | 10 | ho | +1.346361e−02 | 975.27351 | 1.359628e−12 | +5.037384e−11 |
| 0.30 | 10 | vmc-ft | +1.562193e−02 | 323.65732 | 4.096953e−12 | +6.922334e−11 |

Edge ratios (Q²-independent for x ≤ 0.03, where the elastic tail's η range is
already saturated):

| x | σ^el_U(vmc)/σ^el_U(ho) | σ^el_T(vmc)/σ^el_T(ho) |
|---|---|---|
| 0.01 | 1.0106 | +1.0852 |
| 0.03 | 1.0220 | +1.1171 |
| 0.10 | **1.1043** | **−0.2882** (Q² = 5; −0.2957 at 2, −0.2860 at 10) |
| 0.30 | 3.0323 | +3.5501 |

`σ^q_U` — the quasi-elastic tail — is **bit-identical** between the edges at
every point (checked to 1e−13): this knob is the elastic C0 shape and nothing
else. `σ^q_U/σ^el_U` moves only because its denominator does.

**These ratios are smaller than the ones the scoping probe reported**
(×1.015 / ×1.030 / ×1.128 / ×3.28 on σ^el_U and ×1.104 / ×1.141 / ×(−0.43) /
×3.85 on σ^el_T). Two differences, both deliberate: the shipped shape is
r-rescaled to ⟨r²⟩_point = 6.0788 fm² (the probe used the raw transform, whose
own ⟨r²⟩ is 0.9 % low, so its F sits ~5 % high at q = 2 fm⁻¹ and ~8 % high at
2.5), and the shipped shape is continued above 3 fm⁻¹ instead of summing the
noisy table out to q = 278 fm⁻¹. The conclusion is unchanged and its size at
x = 0.10 is if anything **understated** by the probe: the sign still flips.

## B1.6 The headline, republished honestly

The old headline was: **"the tensor fraction of the elastic tail changes SIGN
with x"** (`OPEN_ITEMS_SOLUTIONS.md` §9, `phase_C_numbers.md` §8.1c item 1),
with +1.5595e−04 at x = 0.10 quoted as a result.

What the band says:

* **The sign change survives.** Both edges are negative at x = 0.01 and 0.03 and
  positive at x = 0.30. The statement "the tensor fraction of the elastic tail
  changes sign with x" is **not** an artefact of the unfitted shape.
* **Where it changes is not determined.** On `ho` the crossing is between
  x = 0.03 and 0.10; on `vmc-ft` it is between x = 0.10 and 0.30.
* **At x = 0.10 the sign is INDETERMINATE.** +1.5595e−04 on one edge and
  −4.0702e−05 on the other, at every Q² in the table. The published
  +1.5595e−04 may not be quoted as a result; it is a band edge.
* The `fq_scale` band does not rescue it and makes it worse — at x = 0.10,
  Q² = 5 the two knobs together span

  | `fq_scale` | `ho` | `vmc-ft` |
  |---|---|---|
  | 0 | −9.017624e−05 | −8.166218e−05 |
  | 1 | +1.559536e−04 | −4.070157e−05 |
  | 2 | +4.036705e−04 | +3.523912e−06 |

  i.e. **−9.0e−05 … +4.0e−04**, spanning zero on both edges. (σ^el_T is
  quadratic in F_q, which is why this must be RUN and never rescaled — T12.)
* At x = 0.01 and 0.30 the *magnitude* is a band too, ×1.07 and ×1.17 wide, and
  design §2.1's "O(10⁻²) per unit Q_N" is met at x = 0.30 on both edges and
  missed by 20–60× at x ≤ 0.1 on both. That part of the old conclusion stands.
* `ΔA_zz` from the whole tail moves by **−8.2 %** at x = 0.01, Q² = 5
  (+3.5168e−07 → +3.2291e−07), by **+0.6 %** at x = 0.10 (the quasi-elastic
  dilution dominates there) and by **+37 %** at x = 0.30. All of it is still
  three orders below the band half-width 4.4e−04, so §9's third headline —
  *the tail is not the leading RC systematic at the EIC* — is untouched.

## B1.7 The two policy overrides, made explicitly

**(a) `li6.density` was declared "not an input".**
`run_2026-09-02/design_G_cluster_config.md` §3.2 says it is "**not** an input —
it is the independent validation target" for the assembled α–d one-body density
(T6, `tests/test_cluster_config.cpp`). Reading it in `rc` makes it an input, so
that independence is now **partial**, and the scope of the override is exactly
this: `src/core/cluster_config.cpp` still opens the file **nowhere**, so T6's
4 % radius discrepancy remains an independent test of the cluster model; but any
statement that compares *this form factor* with a cluster-model radius is now
circular and may not be made. Recorded in `include/lipolgen/rc.hpp` beside the
constants, not left implicit.

**(b) The same file is called too wrong to set the tensor normalisation and
right enough to set the monopole shape.** The tension is real; the answer is that
these are different multipoles. `docs/OPEN_ITEMS_SOLUTIONS.md` §5 forbids
deriving a ⁶Li **tensor / b₁** input from these wave functions, and the quantity
that fails is the **quadrupole** — WS98's Q(⁶Li) = −0.23(9) fm² against the
measured −0.0818 fm², a factor 3, which is why `LI6_QUADRUPOLE_FM2` and
`LI6_MU_N` are measured moments. What is read here is the **monopole**, the
ℓ = 0 part of the same density, whose own printed normalisation (2.99922 against
Z = 3, 0.03 % low) and second moment (0.9 % low, and rescaled away) are right to
a fraction of a percent. The quadrupole **normalisation** stays measured on both
edges — `C0Shape` moves F_q's *shape* and never its q → 0 limit, which is why
F_q(0) = −65.91414947 is identical on the two edges (B1.2) and T11 passes
untouched.

## B1.8 Cost

| quantity | measured |
|---|---|
| `Pipeline` setup, `--rc off` | 0.0351 s (median of 7) |
| `Pipeline` setup, `--rc tensor-band --rc-c0-shape ho` | 0.1461 s ⇒ `RcModel` **0.111 s** (the published 0.106 s, same machine, re-measured) |
| `Pipeline` setup, `--rc tensor-band --rc-c0-shape vmc-ft` | **1.3279 s** ⇒ `RcModel` **1.293 s**, i.e. **×11.6** |

The cost is the 99-term j₀ sum, evaluated directly at every η node below the
q = 3 fm⁻¹ cut (≈ 47 % of the log-spaced nodes) rather than interpolated. It is
paid once per run and only on the opt-in edge, so it is published rather than
optimised away — an interpolation table would be a third model choice on top of
the two this task already has to argue.

## B1.9 The ⁶Li form-factor dip: two sites, one provenance

Split off as its own item, because it is a contradiction inside the repository
and not a consequence of anything above.

* `tests/test_rc.cpp` T11 gates q₀ ∈ [2.9, 3.3] fm⁻¹ ⇒ \|t\| ∈ [0.328, 0.424] GeV².
* `docs/PHYSICS_CHANNELS.md`'s coherent row stated, **uncited**, that "the real
  ⁶Li form factor dips near \|t\| ≈ 0.31 GeV²" — that is q = 2.82 fm⁻¹, 9.6 %
  below the shipped q₀ = 3.0999 and **outside** the window T11 gates.

**Determination.** The 0.31 GeV² is the one withdrawn, for three reasons.
(i) It has no source anywhere in the repository — the reference column of that
row is "—" — while the window is at least traceable to
`design_C_tensor_rc.md` §2.1's own bracket around its own admitted starting
guess. (ii) The only ab-initio evidence this repository owns disfavours 2.82
**more** than it disfavours 3.10: the j₀ transform of `li6.density` is
2.454e−03 ± 1.76e−04 at q = 2.82 fm⁻¹ (**14σ** from zero) against
1.076e−03 ± 1.60e−04 at 3.0999 (**6.7σ**), and it has no zero below q ≈ 4.3 fm⁻¹
at all — where its own error bar overtakes it. Even a fully correlated ±1σ shift
of every density bin moves F(3.0999) only to [8.94e−04, 1.26e−03]. (iii) The
argument that sentence exists to make is unaffected either way:
`COHERENT_T_MAX_DEFAULT` = 0.2 GeV² sits below the **whole** gated window
(0.2 < 0.328), which T11 now asserts. *Note 2026-09-23: reason (i) is superseded — 0.31 GeV² ⇔ q = 2.82 fm⁻¹ is Li, Sick, Whitney, Yearian, NPA 162 (1971) 583's q² = 8 fm⁻² minimum (`../../benchmarking/03_data_nuclear.md` §2), and the UVa FB density puts the C0 zero at 2.694 fm⁻¹; see `../run_2026-09-23/phase_B2_nuclear.md` §3. Reasons (ii) and (iii) and the withdrawal as a record of its date stand.*

**Both sites now carry one provenance** — `rc.hpp`'s `LI6_FF_HO_A_FM` /
`LI6_FF_HO_ALPHA` block and T11 — and both now say the same thing, including
that it is a **model** number on either side: q₀ = 3.0999 fm⁻¹ ⇒
\|t\| = (q₀ħc)² = **0.3742 GeV²**, gated over [2.9, 3.3] fm⁻¹, from a starting
guess that no data in this repository can refit, and contradicted at 6.7σ by the
repository's own ab-initio density. The `vmc-ft` edge is the honest expression of
that last fact: it has no zero anywhere.

**Not claimed.** That the real ⁶Li form factor has its first minimum at
3.0999 fm⁻¹, or at 4.3, or anywhere else. No elastic data are in this repository.
What is claimed is that the two sites no longer contradict each other and that
the number they share is labelled as the model number it is.

## B1.10 How to reproduce

```bash
cd /home/cpeng/Projects/polli/LiPolGen && source env.sh
cmake --build build -j
build/lipolgen_tests                       # T11 carries both edges
python -m pytest python/tests/test_rc.py -q
```

Every σ row above, from the shipped Python API alone:

```python
import lipolgen as lg
for edge in ("ho", "vmc-ft"):
    cfg = lg.make_config(isotope="6Li", channel="inclusive", config=1,
                         events=2000, seed=99, rc="tensor-band",
                         rc_c0_shape=edge)
    m = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6)).rc_model
    for x in (0.01, 0.03, 0.10, 0.30):
        u, t, qe = m.tail_sigma_at(x, 5.0)         # per nucleon, GeV^-2
        r0 = m.tail_ratio_at(x, 5.0, 0.0)          # r_U
        r1 = m.tail_ratio_at(x, 5.0, 1.0)          # r_U + r_T
        print(edge, x, (t / 6) / u, qe / u, r0, r1 - r0)
```

`ΔA_zz = (A_zz + 2 r_T)/(1 + r_U) − A_zz` at §8.1c's own `A_zz(Born)` column.
The form-factor rows come from `HoSpin1FF.for_ion(lg.li6(), o)` with
`o.fold_nucleon = False` and `o.c0_shape` set, at `t = (q·HBARC_GEV_FM)²`.
The `F ± σ` table of B1.4 is the transform and its `DRHORP` propagation, done
directly on `data/vmc/density/li6.density`; the continuation variants of B1.4
were run from a scratch `Spin1ElasticFF` subclass that re-implements the
transform independently and reproduces the shipped edge to every printed digit.

## B1.11 What Phase B1 leaves open

* **Q1 is still open.** No fit, no χ²/ndf, no elastic data. The band is a price
  tag, not a measurement.
* **Q10 is still open.** `F_mag`'s (q_z, b) are untouched by this task, and
  `tail_tensor_scale` still bands only their normalisation.
* **The `vmc-ft` edge is not "the right answer".** It is an AV18+UX VMC
  point-proton density whose own ⟨r²⟩ is 0.9 % below the measured one (rescaled
  away here), whose quadrupole is 3× the measured one (not used here), and whose
  transform is noise above q ≈ 4.3 fm⁻¹ (continued here). It is a *second*
  defensible shape, which is what a band needs.
* **Bin-to-bin correlations in the VMC density are not in the ±σ of B1.4.** The
  propagation treats the `DRHORP` column as independent per bin, which is what
  the file supports; a correlated walk would widen it, and would widen the
  "no zero below 4.3 fm⁻¹" statement's error bar with it — never narrow it.

---

# B2 — the s- and p-peaks ship, and `rc_tail` becomes a band of two models

## B2.0 What this task did, in two steps

**Step 1 — promotion, no physics change.** Six functions that lived in an
anonymous namespace in `tests/test_rc.cpp` moved into `src/core/rc.cpp`,
declared in `include/lipolgen/rc.hpp` immediately after `polrad_sigma_qe_u`:

| function | what it is |
|---|---|
| `ll_radiator(z, q2)` | the Weizsäcker–Williams radiator `D(z) = (α/π) ln(Q²/m_e²)(1+z²)/(1−z)` |
| `rosenbluth_spin1(ff, q2, m)` | POLRAD Eq. (A.4) at `Q_N = 0`, spin 1 → `RosenbluthAB{a, b}` |
| `rosenbluth_nucleon(proton, q2)` | the same for a nucleon, through the shipped `nucleon_ff` |
| `dsigma_el_dq2(ab, q2p, sp, m)` | `dσ_el/dQ′²` from `(A, B)` |
| `ll_peaks_spin1(ff, x_a, y, s_a, m)` | the elastic s- and p-peaks → `LlPeaks{s, p}` |
| `ll_peaks_qe(z, n, x, y, s, kf_gev)` | the quasi-elastic ones, `Z` protons + `N` neutrons |

Three deliberate departures from the moved code, each one argued rather than
silent:

1. **Named structs, not `std::pair` and not out-parameters.** `rc.hpp` includes
   `<array> <cstddef> <functional> <memory> <string> <vector>` and **not**
   `<utility>`, and `TailTriple` / `NucleonFF` / `EtaLimits` are how every other
   tuple in that header is carried. `rosenbluth_*` also stopped writing through
   `double*` out-parameters, which have no place in a public header.
2. **Six names are declared, not two.** The task asked for the two entry points.
   T8(c)'s *first* subcase gates `dsigma_el_dq2` at 1e−10 against the laboratory
   Rosenbluth built from `σ_Mott` by a different algebraic route, and it calls
   `rosenbluth_nucleon` to do it. Declaring only the entry points would have
   forced either deleting that gate or keeping a second copy of those two
   functions in the test — and `docs/CONVENTIONS.md` forbids the second copy.
   `ll_radiator` and `rosenbluth_spin1` are declared for symmetry with their
   pairs.
3. **`ll_radiator` returns 0 at `Q² ≤ m_e²`.** Below the electron mass
   `ln(Q²/m_e²)` is **negative**, i.e. a negative cross section, not a small
   one. **This guard is defensive and does not bite in any shipped
   configuration** — measured, the lowest tail-table node is `Q² = 1.5920e−03
   GeV²` for ⁶Li config 1 and `3.9600e−03` for config 2, i.e. **6097×** and
   **1.5e4×** above `m_e² = 2.611199e−07`. (An earlier draft of this section
   claimed the branch was "live in production"; it is not, and the claim is
   withdrawn.) What *could* reach it is a **user** scenario:
   `build_tail_tables` floors `y` at `max(1e−6, y_min)`, so at the sampler's
   `x_min = 1e−4` and config 1's `s = 3980` any `y_min` below **6.6e−4**
   crosses over — the generator scenario's `0.004` clears it by 6.1. T8(c)
   never probes it either, so no gated number moves.

`ll_peaks_qe` also gained a `kf_gev` parameter **defaulting to 0**, the same
signature and the same default as `polrad_sigma_qe_u`. Reason: the two land in
**one** numerator in `tail_ratio_at`, and suppressing the t-peak's
quasi-elastic piece while leaving the s-/p-peaks bare would be an inconsistency
inside a single sum. It is not an extension of POLRAD Eq. (44) but the same
`pauli_suppression(q, k_F)` at each peak's own elastic-vertex `q`. Size: **1.000
at every ⁶Li EIC node with `Q² ≥ 20 GeV²`** (the s-peak's own `Q′²` is far
above `2k_F` there) and **0.89609** at `(x, y) = (0.001, 0.985)`, where
`Q′² = 0.05886 GeV²` gives `q = 0.24464 GeV` against `2k_F = 0.338`. The default 0 keeps every T8(c) number bit for bit.

**Step 2 — the opt-in model.** `RcTailModel::TPeakPlusLL` (`= 2`), reachable as
`--rc-tail-model t-peak+ll`, `make_config(rc_tail_model=…)`,
`lg.RC_TAIL_MODELS`, and `RcOptions.tail_model`. `TailTriple` grew from
`{u, t, qe}` to `{u, t, qe, u_sp, qe_sp}`; `build_tail_tables` tabulates the two
new columns and `tail_ratio_at` adds them to the numerator.

## B2.1 The physics trap, and what was done about it

**There is no tensor s/p peak.** `ll_peaks_spin1` returns the *unpolarised*
Rosenbluth `(A, B)` only; POLRAD supplies no `σ_T` counterpart at the s- or
p-peak. `tail_ratio_at` builds the tensor term as `(q_n/6)·ratio_t·su` with
`ratio_t = σ_t/σ_u`, so folding the s+p into `su` before that product would have
made the tensor tail grow with them — silently asserting that the s-/p-peaks
carry the t-peak's tensor-to-unpolarised ratio. That is an uncited claim, and
`docs/CONVENTIONS.md`'s "no second definition of a physics number" rule makes
deriving one worse than leaving it out.

**What was done:** `u_sp` and `qe_sp` are their own columns and enter the
numerator **outside** the tensor term:

```
num = su + (q_n/6)·ratio_t·su + su_sp + qe_suppression·(sqe + sqe_sp)
                └── unchanged, bit for bit ──┘
```

**The consequence, stated because it is a physics change and not bookkeeping:**
the denominator grows and the tensor numerator does not, so **`TPeakPlusLL`
lowers the tensor FRACTION of the tail** wherever the s-/p-peaks matter. That is
the honest reading — the tensor part of those peaks is **unknown, not zero** —
and it is why `TPeakPlusLL` is a systematic to run *beside* `TPeak` and never a
replacement. Measured, `(w_tail(Q_N) − w_tail(0))/w_tail(0)` at `Q_N = 600`:

> **CORRECTION, 2026-09-06 (tasks B1 + B2).** Two sentences above are now
> wrong. (1) *"POLRAD supplies no σ_T counterpart at the s- or p-peak"* is true
> of **Eq. (38)** and **false of the paper**: Eq. (18) + Eq. (A.4) carries the
> s-/p-peaks' own tensor content and `RcTailModel::PolradFull`
> (`--rc-tail-model polrad-full`, shipped 2026-09-06) **computes** it. (2) the
> tensor part of those peaks is no longer *"unknown, not zero"*: on
> `TPeakPlusLL` it is **bounded, not computed**, by `RcOptions::sp_tensor_scale`
> (B1), and on `PolradFull` it is **computed** — which is why `sp_tensor_scale`
> is refused there for the **opposite** reason it is refused on `TPeak`: the
> term **ran**. Everything else in this section — the `u_sp`/`qe_sp` columns,
> their position outside the tensor term, and every number in the tables — is
> unchanged and still describes `TPeakPlusLL`. See
> `../run_2026-09-06/phase_B_numbers.md` §B1 and §B2.

| x | Q² | `t-peak` | `t-peak+ll` | factor |
|---|---|---|---|---|
| 0.01 | 3 | −0.238354 | −0.0339332 | 0.142 |
| 0.01 | 8 | −0.238354 | −0.226414 | 0.950 |
| 0.10 | 3 | +0.0249065 | +6.45183e−06 | 2.59e−4 |
| 0.10 | 8 | +0.0250016 | +0.00090461 | 0.036 |

(from `build/lipolgen_tests -tc="T8(d)*" -s`, subcase (d); T8(d)(c) gates the
tensor *numerator* itself as unchanged, at 1e−9 — a tolerance set by the
cancellation, see B2.6.)

**ADDED 2026-09-04 — the worst factor is not in this table and is not at low
`Q²`.** At the `Q² ≥ 20 GeV²` cell §B2.3(i) now identifies (`x = 0.7943`,
`y = 0.0088`, `Q² = 27.81`) the same quantity reads **−1.41817e−08 →
−2.20067e−12, factor 1.55e−04** — below every entry above, and *inside* the
window this document previously called safe to 0.16 %. The four rows above are
all at `Q² ∈ {3, 8}`, i.e. outside it. The **absolute** tensor term is
unchanged there as everywhere (both models give −6.975e−16); it is the
denominator that grew by ×6444.

**A tensor s/p peak is not a future refinement of this model.** It needs
Eq. (A.4) evaluated at the shifted `Q′²`, which is `RcTailModel::PolradFull`'s
job and is still **not implemented** and still **throws**.

## B2.2 Mixed approximation orders — stated, not hidden

`TPeak` is POLRAD's η_A quadrature of Eqs. (37)–(39), (43): exact in the peaking
approximation's t channel, integrated over the whole photon phase space it
covers. The s-/p-peaks are a **single-z collinear leading log**: one `D(z)` at
the one `z` elasticity fixes, with no soft exponent, no non-log O(α) term and no
second emission.

**Their sum double-counts nothing** — the three peaks are disjoint regions of
the photon angle — **but it is accurate to the worse of the two**, i.e. to the
leading log, `~1/ln(Q²/m_e²) ~ 5–10 %`. It is a **STATED MODEL, not a
controlled expansion**, and nothing in this repository claims O(α) completeness
for it.

**And its soft `1/(1−z)` is uncancelled.** Both peaks sit at `z → 1` as
`y → 0` (`z_s = (1−y)/(1−x_A y)`, `z_p = 1 − y(1−x_A)`), where `D(z)` diverges.
At `x = 0.744`, `y = 0.0071`, `Q² = 21.04` the radiator is `D(z_s) = 13.50`
(`z_s = 0.993775`, `L = ln(Q²/m_e²) = 18.205`) — far outside the regime where a
single emission is the right expansion. The absolute
contribution there is small (`w_tail − 1 = 3.7e−4`, against `4.5e−8` for the
t-peak alone), so it does not distort the run; but the model is **not
trustworthy in that corner** and no cutoff was invented to hide it.

> **CORRECTION, 2026-09-04 — this paragraph refutes §B2.3(i) below and the
> 2026-09-03 run did not notice.** The point just quoted has `Q² = 21.04 ≥ 20`
> and `y = 0.0071 ≤ 0.9`, and `3.7e−4 / 4.5e−8` is a factor **8200**. It sits
> *inside* the window §B2.3(i) then concluded the two models "agree to < 1 %".
> Reproduced independently: `4.5009e−08 → 3.5749e−04`, **×7942.5**. §B2.3(i) is
> rewritten below; the acceptance survives **only** as an event-weighted
> statement, and the ratio quoted in this paragraph is the FIRST measured
> counterexample to the per-cell one. A `1/(1−z)` that "does not distort the
> run" is a statement about the run's *mean*, never about a cell.

## B2.3 Acceptance, measured — including the criterion that is NOT met

Command: `$SCRATCH/final_b2`, a standalone probe linked against
`build/libLiPolGenCore.so` (source in B2.7).

**(i) ⁶Li EIC config 1, `Q² ≥ 20 GeV²`.** The acceptance asked for agreement to
**< 1 %**.

> **REWRITTEN 2026-09-04.** The table this section carried on 2026-09-03
> sampled `y ∈ {0.5, 0.7, 0.9, 0.985}` and **never probed `y → 0`**. Its
> conclusion — "met at `y ≤ 0.7` (worst 0.16 %), not met above it" — was
> therefore a statement about four hand-picked rows read as a statement about a
> window, and §B2.2 above already contained its counterexample. The rows below
> are unchanged; the low-`y` block and the census are new, and the conclusion
> is replaced.

The original `y ≥ 0.5` rows, reproduced:

| x | y | Q² | t-peak | t + s + p | ratio |
|---|---|---|---|---|---|
| 0.030 | 0.500 | 59.70 | 2.170393e−06 | 2.170513e−06 | 1.00006 |
| 0.100 | 0.500 | 199.00 | 6.237673e−08 | 6.237682e−08 | 1.00000 |
| 0.300 | 0.500 | 597.00 | 1.644884e−09 | 1.644884e−09 | 1.00000 |
| 0.010 | 0.700 | 27.86 | 5.129536e−05 | 5.137751e−05 | **1.00160** |
| 0.010 | 0.900 | 35.82 | 1.425913e−04 | 1.446765e−04 | **1.01462** |
| 0.030 | 0.900 | 107.46 | 8.768410e−06 | 8.774661e−06 | 1.00071 |
| 0.100 | 0.900 | 358.20 | 2.520023e−07 | 2.520073e−07 | 1.00002 |
| 0.010 | 0.985 | 39.20 | 9.414110e−04 | 1.498681e−03 | **1.59195** |
| 0.100 | 0.985 | 392.03 | 1.663755e−06 | 1.682484e−06 | 1.01126 |

**… and the `y → 0` rows the 2026-09-03 table was missing.** All nine sit on
`Q² = 23.88 GeV²` (`x·y = 6.0e−3`), i.e. *inside* `Q² ≥ 20`, walked down in `y`.
`D(z_s)` is the leading-log radiator at the s-peak's own `z`, and `t_min =
M_A²x_A²/(1−x_A)` is where the t-peak's elastic vertex starts:

| x | y | t-peak | t + s + p | ratio | z_s | D(z_s) | t_min [GeV²] |
|---|---|---|---|---|---|---|---|
| 0.01 | 0.6000 | 4.094215e−05 | 4.099824e−05 | 1.0014 | 0.4004 | 0.082 | 8.7e−05 |
| 0.03 | 0.2000 | 1.779724e−06 | 1.783457e−06 | 1.0021 | 0.8008 | 0.351 | 7.9e−04 |
| 0.05 | 0.1200 | 4.193890e−07 | 4.224537e−07 | 1.0073 | 0.8809 | 0.635 | 2.2e−03 |
| 0.10 | 0.0600 | 4.999731e−08 | 5.303204e−08 | **1.0607** | 0.9409 | **1.359** | 8.9e−03 |
| 0.20 | 0.0300 | 5.369098e−09 | 8.789409e−09 | **1.637** | 0.9710 | 2.850 | 0.036 |
| 0.30 | 0.0200 | 1.320914e−09 | 5.251341e−09 | **3.976** | 0.9810 | 4.393 | 0.083 |
| 0.50 | 0.0120 | 1.167129e−10 | 5.652918e−09 | **48.43** | 0.9890 | 7.649 | 0.238 |
| 0.72 | 0.0083 | 4.421547e−12 | 9.923049e−09 | **2244** | 0.9927 | 11.52 | 0.513 |
| 0.80 | 0.0075 | 8.891194e−13 | 1.389794e−08 | **15631** | 0.9935 | 13.00 | 0.644 |

**The census over the sampler's own accepted cells** (the 101 × 77 CLI grid,
3051 accepted cells; `RcModel::tail_ratio(cell, 0)` on both models, i.e. the
ratio the shipped weight actually carries):

| window | cells | > 1 % | by σ | worst | where |
|---|---|---|---|---|---|
| `Q² ≥ 20`, `y ≤ 0.9` | 1356 | **331 (24.4 %)** | **28.2 %** | **×6444** | x = 0.7943, y = 0.0088, Q² = 27.81 |
| … of those, at `y < 0.1` | 440 | 318 (72.3 %) | 82.8 % | ×6444 | same cell |
| `Q² ≥ 20`, `0.15 ≤ y ≤ 0.7` | 660 | **0** | 0 | ×1.0055 | x = 0.0316, y = 0.161 |
| `Q² ≥ 20`, `0.15 ≤ y ≤ 0.8` | 725 | **0** | 0 | ×1.0077 | x = 0.0066, y = 0.771 |
| `Q² ≥ 20`, `0.15 ≤ y ≤ 0.9` | 781 | 5 (0.6 %) | – | ×1.0197 | x = 0.0060, y = 0.845 |

**THE CONCLUSION, in two statements that must never be quoted for each other:**

* **EVENT-WEIGHTED — the acceptance is MET.** Over `Q² ≥ 20`, `y ≤ 0.9`
  (5182 of 200 000 events, seed 1234, P_z = 0 — the suites' fill) the mean dilution `⟨w_tail − 1⟩` moves
  **8.48914e−03 → 8.54072e−03**, **+0.61 %**. Cross-section-weighted over the
  accepted cells: 8.31257e−03 → 8.35972e−03, **+0.57 %**. A rate analysis in
  that window is unaffected at the 1 % level.
* **PER CELL — the acceptance is NOT met, and the failure is at `y → 0`.**
  24.4 % of the window's cells (28.2 % of its cross section) move by more than
  1 %, the worst by **×6444**, and 318 of the 331 sit at `y < 0.1`. The
  per-cell agreement statement that survives is `Q² ≥ 20 GeV²` **and**
  `0.15 ≤ y ≤ 0.7`, **worst 0.55 %** — and even that is *worse* than the
  "0.16 %" the 2026-09-03 table claimed for `y ≤ 0.7`, because that table's
  four `y ≤ 0.7` rows were all at `y ≥ 0.5`.

**Why — and it is NOT the `y → 1` mechanism.** At `y → 1` the s-peak beats
`Y₊ = [1+(1−y)²]/(1−y)` because `z_s = (1−y)/(1−x_A y) → 0` drags the elastic
vertex to `Q′² = z_s Q² → 0`, where the form factor is 1. At `y → 0` the
opposite happens: `z_s → 1`, `Q′² → Q²`, and ⁶Li's coherent form factor is
**dead** there — `ll_peaks_spin1` returns **exactly 0.0** at the worst cell, so
the entire excess is *quasi*-elastic. What wins is the **uncancelled soft
radiator**:

```
1 − z_s = y(1 − x_A)/(1 − x_A y),   1 − z_p = y(1 − x_A)
D(z)    = (α/π) ln(Q²/m_e²) (1+z²)/(1−z)  →  (2α/π) ln(Q²/m_e²) / [y(1 − x_A)]
```

i.e. `D ∝ 1/y`, while the t-peak gets **no** matching growth (`Y₊ → 2` as
`y → 0`) and is simultaneously crushed: the `Q² ≥ 20` cut forces
`x·y ≥ 5.0e−3`, so **low `y` means high `x`**, and the t-peak's own elastic
vertex starts at `t_min ∝ x_A²` — 0.63 GeV² at `x = 0.79`, where the nucleon
dipole `G_D² = (1+t/0.71)^−4` is already down to 0.079 and still falling like
`1/t⁴` — while the s-/p-peaks' vertex stays at `Q′² ≈ Q²`, fixed by the
bin. The table above shows the two effects resolved: **the ratio leaves 1 %
exactly where `D(z_s)` passes 1**, between `x = 0.05` (D = 0.63, ratio 1.0073)
and `x = 0.10` (D = 1.36, ratio 1.061).

**So this corner is a BREAKDOWN OF THE UPPER EDGE, not evidence that `TPeak`
is low by ×6444.** `D(z_s) > 1` means the expansion parameter of a single-`z`
leading log has exceeded one; no single-emission model has an accuracy
statement there. Both tails are minute — at the worst cell `w_tail − 1` is
4.92e−08 (`TPeak`) and 3.17e−04 (`TPeakPlusLL`) — which is exactly why the
event-weighted mean survives. What does **not** survive is the **tensor
fraction of the tail**: `−1.41817e−08 → −2.20067e−12`, **×1.55e−04**, below the
smallest entry of §B2.1's table (2.59e−04, and that one is at `Q² = 3`) and
**inside** the `Q² ≥ 20` window this section used to declare safe. The
*absolute* tensor term is unchanged (both models give −6.975e−16); only the
denominator moved.

**T8(d)(i)** now gates all of it — the raw-quadrature ratio at the worst point,
the cell census, the `D(z_s) > 1 ⟺ ratio > 1 %` correspondence on the
`Q² = 23.88` line, the `0.15 ≤ y ≤ 0.7` window and the tensor collapse — so the
corner cannot regress into invisibility again. T8(d)(f) still asserts the
`y = 0.985` failure. On the test binary's coarser 40 × 24 grid the census reads
**45 of 188 cells (23.9 %), 44 of them at `y < 0.1`, worst ×5667** — the same
statement at a different resolution.

**(ii) The HERMES deuteron point — met, and it reproduces T8(c) exactly.**
`x = 0.012`, `y = 0.85`, `Q² = 0.5283`, T8(c)'s own unsuppressed convention:

| quantity | value |
|---|---|
| t-peak alone | 3.603076e−03 |
| t + s + p | 1.571906e−02 |
| **ratio** | **4.3627** (acceptance: ≥ 4×) |
| t-peak fraction | 0.229217 — the same six digits T8(c) printed before the move |
| … elastic / quasi-elastic | 0.649504 / 0.165690 |

**(iii) The low-Q² corner — met, and it corrects a documented wording.**
`x = 0.01`, `y = 0.1`, `Q² = 3.98`, quasi-elastic:

| `k_F` | t-peak | s-peak | p-peak | s/t | (s+p)/t |
|---|---|---|---|---|---|
| 0 (T8(c)'s convention) | 1.324175e−05 | 3.090839e−05 | 1.341732e−05 | **2.3342** | **3.3474** |
| 0.169 (the shipped default) | 6.251978e−06 | 3.090839e−05 | 1.341732e−05 | **4.9438** | **7.0899** |

The documented **3.3×** is reproduced to three digits — but it is **(s+p)/t**,
not "the s-peak", which is 2.33×. Four documents said "the quasi-elastic
**s-peak** is already 3.3× the t-peak"; all four now say s+p, give the s-peak
alone, and give the shipped-`k_F` numbers beside the `k_F = 0` ones.

**(iv) `--rc off` is byte-identical.** And so is the default `--rc
tensor-band`:

```
$ git archive ac22331 | tar -x -C $SCRATCH/base_tree     # untouched baseline
$ cmake -S $SCRATCH/base_tree -B $SCRATCH/base_tree/build -DCMAKE_BUILD_TYPE=Release && cmake --build … -j
$ for t in "" base; do python -m lipolgen.cli --isotope 6Li --channel inclusive \
      --config 1 --events 5000 --seed 1234 --rc off --npz …; done
```

| comparison | arrays | result |
|---|---|---|
| `--rc off`, `ac22331` vs this tree | 48 | **byte-identical**, `meta` identical |
| `--rc tensor-band` (default `t-peak`), same | 51 | **byte-identical**; `meta` differs only by B1's `rc_c0_shape` key and its provenance sentence |

The comparison is on the raw `float64` bytes (`x.view(np.uint8)`), not a
tolerance. The structural reason it holds: under `TPeak` the two new tables are
identically zero and `tail_ratio_at`'s s+p block is **not entered at all**, so
not one floating-point operation of the shipped default moved. T8(d)(a) gates
that from inside the test binary.

## B2.4 The clip statistics — regenerated, and the surprise

The task expected the new tail to move `clipped_cell_fraction()`,
`clipped_fraction_by_y()` **and** the per-event `kRcClipTail` rate. **The first
two move; the third does not, on either grid or at either event count.**

| statistic | `t-peak` | `t-peak+ll` |
|---|---|---|
| `clipped_cell_fraction()`, 101 × 77 CLI grid | 0.000000 | **0.003600** |
| `clipped_fraction_by_y()` — `y < 0.5` / `0.5–0.9` / `> 0.9` | 0 / 0 / 0 | **0 / 0 / 0.046205** |
| `clipped_cell_fraction()`, T8(d)'s 40 × 24 test grid | 0.000000 | **0.010093** |
| `clipped_fraction_by_y()` on that grid | 0 / 0 / 0 | **0 / 0 / 0.058537** |
| `meta["rc_clipped_tail_events"]`, 2 000-event CLI run | 0 | **0** |
| … at 200 000 events | 0 | **0** |

**Why the event count stays 0 while 4.6 % of the `y > 0.9` nodes clip:** the
node statistic is evaluated on the table's own `(x, y)` node grid, whose `y`
axis carries refinement rows at `y ∈ {0.9, 0.95, 0.97, 0.98, 0.985}` and whose
`x` axis is the sampler's cell **edges** — down to `x = 1e−4`. The events land
at accepted cell **centres**, none of which reaches the `(low x, y → 1)` corner
where the ratio exceeds `tail_max = 10`. This is precisely the divergence
`rc.hpp` and `PHYSICS_CHANNELS.md` already warn about ("nodes are not
event-weighted"), now with a run where the two actually disagree.

**Nothing in the shipped default's pinned numbers had to be regenerated.** The
`0 % / 0 / 0 / 0` node fractions and the `0 / 0` event counts in
`phase_C_numbers.md` §8.3, `OPEN_ITEMS_SOLUTIONS.md` and `PHYSICS_CHANNELS.md`
are still exactly right for `TPeak`; what the documents gained is a
`t-peak+ll` column beside them, and the sentence saying which model each number
belongs to.

## B2.5 What the whole run does

200 000 inclusive ⁶Li config-1 events, `seed = 1234`, `tensor-thirds(0, 0.6)`:

| | `t-peak` | `t-peak+ll` |
|---|---|---|
| mean `rc_tail` | 1.02183630 | **1.04036893** |
| median `rc_tail` | 1.000214 | 1.007246 |
| max `rc_tail` | 3.410153 | 5.981595 |

The tail **dilution** (`mean − 1`) roughly doubles, 2.18 % → 4.04 %, driven
almost entirely by the `Q² < 20 GeV²` bulk of the scenario (the run spans
`Q² = 1 … 986`; only 2.7 % of events sit above 20). Kinematics, `weight` and
every non-`rc` column are **bit-identical** between the two runs —
`RcModel::fill` takes no `Rng&`, and the pytest gate asserts
`np.array_equal` on `x`, `q2` and `weight`.

Per-event `(w_tail − 1)` ratios are **not** a useful summary and are recorded
here only so nobody re-derives them and panics: they reach 9.3e4, because at
high `x` and low `y` the t-peak is `~1e−8` of the Born while the s+p is
`~4e−4` — a huge *relative* change of a quantity that is still four orders of
magnitude below anything that matters.

## B2.6 The seven-file chain, and the two tests that had to change

| file | change |
|---|---|
| `include/lipolgen/rc.hpp` | six declarations + `LlPeaks` / `RosenbluthAB`; `RcTailModel::TPeakPlusLL` with the three-paragraph trap comment; `rc_tail_model_name`; `TailTriple` → 5 wide; `sigma_tail_u_sp()` / `sigma_tail_qe_sp()`; the two `LOWER BOUND` blocks rewritten |
| `src/core/rc.cpp` | the six definitions; `rc_tail_model_name`; the constructor guard now refuses `PolradFull` only; `tail_sigma_at`, `build_tail_tables`, `tail_ratio_at` |
| `src/core/pipeline.cpp` | `validate()` accepts `TPeakPlusLL`; the `m_lepton` cross-check no longer keys on `== TPeak` |
| `python/bindings.cpp` | `meta["rc_tail_model"]` now calls `rc_tail_model_name`; the `TPeakPlusLL` enum value; two new table properties; `tail_sigma_at` returns 5; three docstrings |
| `python/lipolgen/__init__.py` | `RC_TAIL_MODELS`, `__all__`, `make_config(rc_tail_model=…)` |
| `python/lipolgen/cli.py` | `--rc-tail-model`, `DEFAULTS`, the `make_config` call, and a **model-aware run banner** (it hard-coded "the tail is the t-PEAK ONLY") |
| `tests/`, `python/tests/` | T8(c) calls the shipped names; **T8(d)** is new (9 subcases, the ninth added by §B2.9); `test_rc_pipeline.cpp`'s "tail model" subcase now also asserts `TPeakPlusLL` validates, plus a new `m_lepton` subcase; `test_rc.py` gains `test_tail_model_knob_is_a_seven_file_chain` |

**`m_lepton` is still refused, and the reason changed.** `TPeak` carries no
lepton mass at all. `TPeakPlusLL` **does** — `ll_radiator` is
`(α/π) ln(Q²/m_e²)(1+z²)/(1−z)` — but it reads `constants.hpp`'s `M_ELECTRON`
directly and **not** `rc_options.m_lepton`, deliberately: a muon radiator would
also need the muon's own elastic kinematics, so honouring `m_lepton` in the log
alone would be a half-change dressed as a whole one. `validate()` therefore
refuses it on **both** implemented models, and the rule stays intact — what
`meta` records is what ran.

**One tolerance is 1e−9 and not 1e−12, and that is not slack.** T8(d)(c)
compares `tail_ratio_at(x, q2, 600) − tail_ratio_at(x, q2, 0)` between the two
models. Adding the s+p makes the *ratio* up to 3900× larger at `x = 0.1,
Q² = 3` while the tensor term does not move, so the subtraction loses ~5 digits
and the double-precision floor is `ε·base/|d| ≈ 2e−16 · 1.6e5 ≈ 3e−11`.
Measured worst case over the four gated points: **1.47e−11**; the other three
are 1.2e−15, 2.0e−14 and exactly 0. A 1e−12 gate there would have been gating
round-off.

**One new gate is not about `TPeakPlusLL` at all.** T8(c) now checks that
`rosenbluth_spin1` is **not** a second definition of Eq. (38)'s unpolarised
integrand: that integrand *is* `A·X̃ − B/(2η)` identically, so running the
shipped `(A, B)` through the shipped `polrad_tpeak_quadrature` must reproduce
the shipped `polrad_sigma_el_u`. Gated at **1e−14** over 9 `(x_A, y)` points.

## B2.7 How to reproduce

```bash
source env.sh
cmake --build build -j
build/lipolgen_tests -tc="T8(c)*" -s          # the promoted construction
build/lipolgen_tests -tc="T8(d)*" -s          # the new model, 9 subcases
python -m pytest python/tests/test_rc.py -q

# The B2.5 / B2.3(i) PAIR.  `--pz 0` is LOAD-BEARING and `--rc-tail-model` has
# to be run on BOTH edges: the published run is `tensor_thirds_plan(0, 0.6)`
# (B2.5), while the CLI's own default is pz = 0.7, which fills the spin
# categories differently and therefore draws a different event sample.
for M in t-peak t-peak+ll; do
  python -m lipolgen.cli --isotope 6Li --channel inclusive --config 1 \
      --events 200000 --seed 1234 --pz 0 --pzz 0.6 \
      --rc tensor-band --rc-tail-model "$M" --npz "run_$M.npz"
done
python - <<'PY'
import numpy as np
r = {}
for m in ("t-peak", "t-peak+ll"):
    d = np.load("run_%s.npz" % m, allow_pickle=True)
    w = (d["q2"] >= 20.0) & (d["y"] <= 0.9)
    r[m] = (int(w.sum()), float((d["rc_tail"][w] - 1.0).mean()),
            float(d["rc_tail"].mean()))
print(r)                       # B2.3(i)'s window pair and B2.5's run means
PY
```

**Why `--pz 0` is in that command.** Until 2026-09-04 this block ran the CLI at
its defaults, and the defaults do **not** reproduce the numbers printed beside
it. Measured both ways:

| | window events | `⟨w_tail − 1⟩`, `t-peak → t-peak+ll` | mean `rc_tail`, `t-peak → t-peak+ll` |
|---|---|---|---|
| **`--pz 0`** — the published run | **5182** | **8.48914e−03 → 8.54072e−03** (+0.6076 %, the "+0.61 %" of B2.3(i)) | **1.02183630 → 1.04036893** |
| CLI default `pz = 0.7` | 5194 | 8.71965e−03 → 8.77375e−03 (+0.6205 %) | 1.02177853 → 1.04027419 |

Same seed, same everything else; only the fill's `P_z` differs, and with it
which events the plan's categories draw. The conclusions are identical on both
— that is the point of showing the second row — but the *digits* B2.3(i) and
B2.5 publish are the first row's, and a recipe has to produce the digits it is
printed beside.

The acceptance tables of B2.3 come from a standalone probe compiled with

```bash
g++ -std=c++17 -O2 -I include -I third_party probe.cpp -o probe \
    -L build -lLiPolGenCore -Wl,-rpath,$PWD/build
```

calling `polrad_sigma_el_u`, `polrad_sigma_qe_u`, `ll_peaks_spin1` and
`ll_peaks_qe` directly — no table, no interpolation — with the deuteron form
factor built from T8(c)'s own `(mu_n, q_fm2, a_fm, alpha, fm_qz_fm, fm_b_fm)`
stand-in, and `HoSpin1FF::for_ion(LI6())` for the ⁶Li rows. The byte-identity
comparison of B2.3(iv) reads both `.npz` files with `allow_pickle=True` and
compares `x.view(np.uint8)`.

## B2.7a Suite tallies

| | before B2 (B1 applied) | after B2 | after the 2026-09-04 correction |
|---|---|---|---|
| `build/lipolgen_tests` | 374 cases / **17 209 856** assertions / 1 skipped / 0 failed | **375** cases / **17 215 875** assertions / 1 skipped / **0 failed** | **377** cases / **17 216 856** assertions / 1 skipped / **0 failed** |
| `python -m pytest python/tests -q` | 263 passed | **264 passed** | **268 passed** |
| `validation/check_physics_channels_links.py` | 1061 refs / 0 broken | **1061 refs / 0 broken** | **1071 refs / 0 broken** |

(The 2026-09-04 column is the whole B1–B6 tree with §B2.9's correction applied.
T8(d)(i) added 230 assertions inside the existing `T8(d)` case, so the case
count moves only by what B3–B6 added after B2 was measured.)

(The task brief quoted the *Phase B start* tallies — 374 / 17 208 032 / 261 /
1049 — which is commit `ac22331` before B1. B1 added 1 824 assertions and 2
pytests and 12 doc references; the "before B2" column above is the tree B2
actually started from, measured.)

`--fix` was run on `check_physics_channels_links.py` because this task moved
several hundred lines in `rc.hpp` and `rc.cpp`. **Its diff was reviewed and 21
of its rewrites were reverted by hand.** The fixer resolves a broken reference
to the *nearest* line containing the symbol name, which for a symbol also
mentioned in prose above its own declaration lands on the **comment**, not the
declaration — e.g. `` `include/lipolgen/rc.hpp:675` `qe_suppression` (as of 66dcda2) `` (the
field) became `:592` (a sentence in `nucleon_ff`'s block). All 21 now point at
the declaration again (`qe_suppression` → `:958`, `RcOptions` → `:919`,
`RcModel` → `:1031`, `RcTailModel` → `:828`, `HoSpin1FF` → `:507`,
`Pipeline::event` → `pipeline.cpp:1379`, and so on). Checked mechanically: of
the 148 references whose line moved, **0** now land on a comment line that
previously landed on code.

## B2.8 What Phase B2 leaves open

* **`RcTailModel::PolradFull` is still not implemented and still throws.** It
  is the only route to a **tensor** s/p peak, and therefore the only route to a
  tail whose tensor fraction is trustworthy at fixed-target kinematics or at
  the `y → 1` edge.
* **The soft `1/(1−z)` is not exponentiated.** `TPeakPlusLL` overshoots as
  `y → 0`, where `D(z) > 1` — which, measured on 2026-09-04, is the whole of
  the `Q² ≥ 20 GeV²` window below `y ≈ 0.06`. No cutoff was invented; the
  corner is documented, gated by T8(d)(i) and small in absolute terms, but the
  upper edge has **no accuracy statement at all** there.
* **The `<1 %` acceptance at `Q² ≥ 20 GeV²` is EVENT-WEIGHTED-ONLY** (+0.61 %
  on `y ≤ 0.9`). **Per cell it fails**: 24.4 % of that window's accepted cells
  move by more than 1 % and the worst by ×6444, all at `y → 0`. The surviving
  per-cell window is `0.15 ≤ y ≤ 0.7` at 0.55 %. §B2.3(i) carries the
  measurement; the 2026-09-03 "0.16 % at `y ≤ 0.7`" is withdrawn. Nothing was
  adjusted to make the original claim true.
* **A tensor fraction of the tail quoted in that corner is meaningless.** It
  falls by ×1.55e−04 there *purely because the unpolarised denominator grew*,
  and the numerator the fall is measured against is the t-peak's own.
* **`Mo–Tsai` is still not obtained**, so the *absolute* normalisation of
  neither edge is validated against an exact tail. The band prices the *spread*
  between two approximations, which is not the same thing as an error bar.
* **No polarised quasi-elastic tail, on any peak.** Unchanged by this task.

## B2.9 The 2026-09-04 correction — a published acceptance that was false

**What was wrong.** Four shipped sites (`include/lipolgen/rc.hpp`,
`docs/USAGE.md`, `docs/PHYSICS_CHANNELS.md`, `docs/OPEN_ITEMS_SOLUTIONS.md`)
and §B2.3(i) above said the two tail models "agree to 0.16 % at `y ≤ 0.7` and
1.5 % at `y = 0.9`" for ⁶Li at `Q² ≥ 20 GeV²`. **A reader takes that as
per-cell safety and it is not true**: 24.4 % of the window's accepted cells
disagree by more than 1 % and the worst by ×6444. The failure is at **`y → 0`**,
which the "and `y ≤ 0.9`" qualifier added on 2026-09-03 does not touch, and
§B2.2 above **already contained a counterexample** (`x = 0.744`, `y = 0.0071`,
`Q² = 21.04`, `3.7e−4` against `4.5e−8` — a factor 8200) fourteen lines above
the table that concluded the criterion was met. Nothing in the code was wrong;
the *claim* was, and its *scope* was understated.

**Why no test caught it.** T8(d)(e) checked `y = 0.5`; T8(c)(i) `y ∈ {0.5,
0.9}`; T8(d)(f) `y = 0.985`. The acceptance table sampled `y ∈ {0.5, 0.7, 0.9,
0.985}`. **Nothing probed `y → 0`** — the region the model's own
"uncancelled soft `1/(1−z)` overshoots as `y → 0`" caveat names.

**What changed.** §B2.1, §B2.2 and §B2.3(i) above; the four shipped sites; the
`Q² ≥ 20` sentences in `phase_C_numbers.md` and `run_2026-09-03/STATUS.md`; two
`python/lipolgen/cli.py` help strings that still carried the un-narrowed
sentence (`--rc`: "validated for ⁶Li at Q² ≥ 20 GeV²"; `--rc-tail-model`:
"Negligible for ⁶Li at EIC Q² ≥ 20 GeV² (< 1 %)"); the run banner; and a new
gate, **T8(d)(i)**, which pins the raw-quadrature ratio at the worst point, the
cell census, the `D(z_s) > 1 ⟺ ratio > 1 %` correspondence on the `Q² = 23.88`
line, the surviving `0.15 ≤ y ≤ 0.7` window and the tensor-fraction collapse.
**No physics number moved**; `--rc off` and the `TPeak` default are untouched
(comments, help text, docs and one new test only).

**And one number was unsourceable.** The run banner said the t-peak is "12 % at
the `y → 1` edge". That sentence entered on 2026-09-03 (it is absent at
`ac22331`), appears in **no other document**, and no measurement reproduces it:
the t-peak fraction on the `y = 0.985` row is **62.8 %** at `x = 0.01`, and its
**minimum over the whole `y ≥ 0.9` part of the tail table's own 101 × 77 node
grid** — which tops out at the scenario's `y_max = 0.985` — is **36.7 %**
(`x = 1e−4`, `y = 0.985`, `Q² = 0.392`). Never 12 %. It is replaced by the
measured **62.8 %**, and the banner now states the per-cell census beside the
event-weighted figure.

| where on the `y = 0.985` edge | t-peak / (t + s + p) |
|---|---|
| x = 1e−4 (`Q² = 0.392`) | **36.73 %** ← the minimum over the table's `y ≥ 0.9` nodes |
| x = 1e−3 (`Q² = 3.92`) | 48.94 % |
| x = 0.01 (`Q² = 39.2`) | **62.82 %** ← the banner's new number |
| x = 0.10 (`Q² = 392`) | 98.89 % |

(Push past the table and it keeps falling — 17.4 % at `x = 0.098`, `y = 0.999`
— but no shipped scenario runs there, and no accepted **cell** reaches the
edge's worst corner either: over the 3051 accepted cells with `y ≥ 0.9` the
worst ratio is only **×1.92** (t-peak 52.2 %, at `x = 8.7e−4`, `y = 0.971`).
That is §B2.4's node-vs-event divergence again, and it is a second reason the
`y → 1` edge was never the corner that mattered for this generator's runs — the
`y → 0` one is.)

**How to reproduce the correction.**

```bash
source env.sh && cmake --build build -j
build/lipolgen_tests -tc="T8(d)*" -sc="*(i)*" -s     # the new gate, with its MESSAGEs
```

The census and the event-weighted means come from the shipped Python surface —
`Pipeline.dis_sampler` for the accepted cells and `RcModel.tail_ratio` /
`tail_sigma_at` on two pipelines that differ only in `rc_tail_model`, plus a
200 000-event `seed = 1234` pair for the event-weighted figure. On the test
binary's coarser 40 × 24 grid the same census reads 45 of 188 cells (23.9 %),
44 at `y < 0.1`, worst ×5667.

---

# B3 — the polarised quasi-elastic tail gets a price tag, not a calculation

## B3.0 What this task did, and what it deliberately did not

`polrad_sigma_qe_u` has no tensor partner, and `TailTriple::qe` entered
`tail_ratio_at` with no `Q_N` coefficient at all: the quasi-elastic tail was
treated as **exactly tensor-blind**. After Q9 closed the Pauli factor and the
per-nucleon fix raised the quasi-elastic term to **22 % / 73 % / 99.9 %** of
`rc_tail` at x = 0.01 / 0.10 / 0.30, that zero was the largest unpriced piece
of the tail.

**It is still not computed, and this task did not compute it.** POLRAD supplies
no tensor partner to Eq. (44); no polarised quasi-elastic radiative-tail
calculation exists for an A = 6 spin-1 nucleus; Zhou *et al.*, PRL **82** (1999)
687 — the citation the neglect rests on — is a **deuteron** measurement and is
**not in this repository**. What shipped is a **stand-in**:

```
  num += qe_suppression * (Q_N/6) * qe_tensor_scale * (σ^el_T/σ^el_U) * σ^q_U
```

one term in `RcModel::tail_ratio_at`, reusing the **same** `ratio_t` the
elastic tensor term uses, so `qe_tensor_scale = 1` is literally the sentence
*"the quasi-elastic tensor fraction is the coherent elastic one"*.
**Default 0.0.**

## B3.1 The argument, stated as an argument and not as a bound

The brief called scale = 1 the bound *"the quasi-elastic tensor fraction is at
most the elastic tensor fraction"*. **It is not a bound**, and the header says
so in those words. What it is, and where it breaks:

1. **It borrows a coherent quantity for an incoherent process.** `σ^el_T` is
   built from F_Q and F_M — rank-2 properties of the **whole** ⁶Li ground-state
   charge and magnetisation distributions. `σ^q_U` is an incoherent sum of
   Z + N **nucleon** elastic tails; its tensor partner would come from the
   alignment of the nucleon **momentum** distribution inside an aligned ⁶Li
   (plus the off-shell and Pauli response), a different object with its own Q²
   shape. No step of any derivation connects the two.
2. **Where it is conservative.** It replaces an exact zero on the dominant
   piece of the tail with a number that is reported, banded, linear, and
   carried in the same (x, y) weighting as everything else in the numerator. It
   cannot be quietly wrong the way a zero can.
3. **Where it fails, and in which direction — measured.** ⁶Li's elastic tensor
   fraction is *anomalously small at low x*, and for a reason the quasi-elastic
   piece has no reason to share:

   | Q² = 5 GeV² | x = 0.01 | x = 0.10 | x = 0.30 |
   |---|---|---|---|
   | ⁶Li `σ^el_T/σ^el_U` (this tree, `ho` edge) | **−3.0566e−03** | **+9.3572e−04** | **+7.9721e−02** |
   | deuteron `σ_q/σ_u` (POLRAD's own, transcription check §8; x = 0.10 column re-driven 2026-09-23, was +0.064) | +0.106 | +0.062 | −0.117 |
   | ratio | **35×** | **67×** (68× until 2026-09-23) | 1.5× |

   Design §2.1 says why the ⁶Li number is small: what sets it is Q_A/Z, and
   Q(⁶Li)/Q(d) = 0.29 against Z = 3. But that near-vanishing quadrupole is a
   cancellation in the **coherent charge distribution**, while the alignment of
   the **nucleon momentum distribution** — what a quasi-elastic tensor response
   rides on — is carried by the deuteron-like pair that holds the whole spin.
   So at x ≲ 0.1 **scale = 1 may understate the omission by one to two orders
   of magnitude**, and a run pricing it should quote **both 1 and O(10²)**.
   (The deuteron *elastic* ratio is itself only a proxy for a deuteron
   *quasi-elastic* one, which nobody has computed either. It bounds nothing; it
   says the right scale is not 1.)
4. **The sign is meaningless.** `σ^el_T/σ^el_U` changes sign with x and its
   sign at x = 0.10 is a C0-shape band edge (§B1). The stand-in inherits that
   sign and there is no argument that the quasi-elastic tail shares it. Read
   the **magnitude**.

## B3.2 ΔA_zz from the stand-in — and a correction to the reference it is measured against

`ΔA_zz = (A_zz + 2 r_T)/(1 + r_U) − A_zz`, the shift a fit that ignored the
tail would make. `r_U` is unchanged by this knob **exactly** (the term carries
`Q_N`), so the whole effect is in `r_T` and

```
  Δ(ΔA_zz) = 2 [r_T(scale) − r_T(0)] / (1 + r_U)
```

is independent of `A_zz`. At Q² = 5 GeV², ⁶Li `--config 1`, defaults, `ho`:

| x | `r_U` | `r_T` at scale 0 | `r_T` at scale 1 | tensor tail × | ΔA_zz(0) | ΔA_zz(1) | **Δ(ΔA_zz)** |
|---|---|---|---|---|---|---|---|
| 0.01 | 5.18222e−04 | −2.05867e−07 | −2.63998e−07 | **×1.2824** | −1.76980e−07 | −2.93181e−07 | **−1.16201e−07** |
| 0.10 | 1.31963e−06 | +5.49173e−11 | +2.05801e−10 | **×3.7475** | +2.12745e−09 | +2.42922e−09 | **+3.01768e−10** |
| 0.30 | 6.86176e−08 | +8.93669e−13 | +8.98872e−10 | **×1005.8** | +5.13016e−12 | +1.80109e−09 | **+1.79596e−09** |

The `×1005.8` at x = 0.30 is not a surprise, it is the point: there the tail is
99.9 % quasi-elastic, so lending that piece the elastic tensor fraction
multiplies the tensor part of the tail by `1 + σ^q_U/σ^el_U`.

**Against the band, and here a published reference number has to be corrected.**
The brief asks for the comparison against a band half-width of **4.4e−04**.
That number is `δ(x)·|A_zz(Born)|` from `phase_C_numbers.md` §8.2, and its
`A_zz(Born)` column **does not reproduce**: it is **3.253983×** this
generator's own ⁶Li A_zz at every x and Q² in the table. The cause is
identified, not guessed — that column was built with the **deuteron** Miller b₁
table (`toy_b1`, `B1Mode::Digitized`) instead of the pipeline's own ⁶Li b₁
(`Li6B1(MillerB1)`), which the shipped `default_inclusive_kernel` has used for
every spin-1 ion since before the RC commit. Reproduced to five digits:

| x, Q² = 5 | published `A_zz(Born)` | `azz(toy_b1 …)` — the deuteron table | **this generator's ⁶Li A_zz** |
|---|---|---|---|
| 0.01 | −1.47350e−03 | **−1.473482e−03** | **−4.528242e−04** |
| 0.10 | −4.97510e−03 | **−4.975091e−03** | **−1.528923e−03** |
| 0.30 | −1.58520e−04 | **−1.585230e−04** | **−4.871661e−05** |

So the band half-widths at Q² = 5 are **1.3585e−04 / 9.6800e−05 / 7.3075e−07**
at x = 0.01 / 0.10 / 0.30, not 4.42e−04 at x = 0.01. (`phase_C_numbers.md`
§8.1c/§8.2's `A_zz`, `τ`, `w_hi − 1` and `ΔA_zz` columns are all affected;
their `σ^el_T`, `σ^q_U` and `r_U` columns are **not** and reproduce here to
every printed digit. Re-publishing those tables is **not** part of B3 and was
not done — the finding is recorded at both sites, with its cause, so the next
task can do it deliberately.)

The `ΔA_zz` column of §8.1c does not merely shrink: `ΔA_zz =
[2 r_T − A_zz r_U]/(1 + r_U)`, and at x = 0.01 the two terms are comparable, so
it **changes sign** — published +3.5167e−07 against **−1.76980e−07** measured
here. At x = 0.10 and 0.30 the sign survives and the magnitude falls
(+6.6751e−09 → **+2.12745e−09**, +1.2665e−11 → **+5.13016e−12**).

**The stand-in against the band, both denominators:**

| x | Δ(ΔA_zz) at scale 1 | band half-width (⁶Li b₁, this tree) | fraction | vs the published 4.42e−04 |
|---|---|---|---|---|
| 0.01 | −1.16201e−07 | 1.3585e−04 | **0.086 %** | 0.026 % |
| 0.10 | +3.01768e−10 | 9.6800e−05 | **0.00031 %** | — |
| 0.30 | +1.79596e−09 | 7.3075e−07 | **0.246 %** | — |

**Reading.** At scale 1 the stand-in is a **sub-percent** entry in the RC
systematic budget everywhere, and the band on δ(x) still dominates it by
400–300 000×. It is *largest, relatively*, at **x = 0.30** — where δ(x) is at
its 0.015 floor and the tail is 99.9 % quasi-elastic — and that is the opposite
end of x from where the tail's absolute size peaks. **But the scale is not
known to be 1**: at the O(10²) of §B3.1 the x ≤ 0.1 entries become
8.6 % and 0.03 % of the band and the x = 0.30 entry exceeds it. That is the
honest statement of what was unpriced.

## B3.3 The three departures from the brief, each argued

1. **The term rides inside `qe_suppression`.** The brief specified
   `+ (q_n/6) * opt_.qe_tensor_scale * ratio_t * sqe` appended to the
   numerator. It ships multiplied by `qe_suppression` as well, because that
   knob is documented and tested as *a flat multiplier on the whole
   quasi-elastic tail*: at `qe_suppression = 0` the brief's form would leave a
   quasi-elastic **tensor** term whose unpolarised parent had been switched
   off — a tail with no parent, and it would have broken the reading of T8''s
   "`qe_suppression = 0` reproduces the elastic-only tail" for any caller who
   set both. A doctest subcase now gates both halves of the new behaviour
   (`qe_suppression = 0` kills the stand-in exactly; `0.5` halves it exactly).
2. **`qe_suppression = 0` is *not* refused with a non-zero scale**, though
   it degenerates the same way `with_qe_tail = false` does. `with_qe_tail` is a
   structural switch and is refused (rule 5 of the brief, the `m_lepton` rule);
   `qe_suppression = 0` is a *numeric edge of a knob that did run*, exactly
   like `qe_kf_gev = 0`, which is not refused either. Refusing one and not the
   other would have been a new, inconsistent rule.
3. **A CLI switch was added** (`--rc-qe-tensor-scale`), not only the API knob,
   because "mirror `qe_suppression` through its plumbing sites" includes
   `cli.py`, and every field `make_config` takes has a switch. That made
   `PHYSICS_CHANNELS.md`'s "Six CLI switches" wrong — it had **already** been
   wrong since `--rc-c0-shape` and `--rc-tail-model` landed — and it now reads
   nine, listed.

**And what the brief said not to do was not done:** there is no
quadratic-vs-linear T12 analogue, because the term is one product and is
**exactly** linear. That is asserted instead: scale 2 and 4 give exactly 2× and
4× the increment (rtol 1e-10) at eight (x, Q²) points, and the increment is
linear in `Q_N` too (the −2 : +1 ratio at 1e-9), which is what keeps T14's
Eq. (43) identity exact.

## B3.4 The default did not move — measured, not argued

The numerator was **restructured**, not merely extended, so bit-identity was
checked empirically rather than reasoned about. The one-line change was reverted
in place, the library rebuilt, and two CLI runs regenerated at the same seed:

```
python -m lipolgen.cli --isotope 6Li --channel inclusive --config 1 \
    --events 2000 --seed 4242 --rc {tensor-band,off} --npz …
```

**Every numeric column and every `meta` entry is byte-identical** between the
pre-change and post-change builds, on both `--rc tensor-band` and `--rc off`
(compared as `ndarray.tobytes()` per column plus a key-by-key `meta` diff).
Whole-run mean `rc_tail` on 200 k events, seed 1234, reproduces §B2's published
default to nine digits:

| `qe_tensor_scale` | mean `rc_tail` | min | max | clipped tail events |
|---|---|---|---|---|
| **0 (default)** | **1.021836305** | 1.000000035 | 3.410153365 | 0 / 200 000 |
| 1 | 1.021836323 | 1.000000027 | 3.410062720 | 0 / 200 000 |
| 100 | 1.021838073 | **0.985845619** | 3.401088866 | 0 / 200 000 |

Two things in that table are worth stating. The **mean barely moves** even at
scale 100, because the stand-in carries `Q_N/6` and the run averages over
m = ±1, 0. And at scale 100 the **minimum drops below 1**: the tensor term is
negative where `σ^el_T` is, so a large enough stand-in makes `rc_tail` a weight
that *removes* events at some (x, Q², m). Node clipping stays 0 % globally and
in all three y bands at every scale tested.

## B3.5 How to reproduce

```bash
source env.sh
cmake --build build -j
build/lipolgen_tests -tc="B3*" -s
python -m pytest python/tests/test_rc.py -q
```

Every σ and ΔA_zz row above, from the shipped Python API alone:

```python
import lipolgen as lg
def mk(**kw):
    cfg = lg.make_config(isotope="6Li", channel="inclusive", config=1,
                         events=2000, seed=99, rc="tensor-band", **kw)
    return lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
p0, p1 = mk(), mk(rc_qe_tensor_scale=1.0)
m0, m1 = p0.rc_model, p1.rc_model
k, s = p0.dis_sampler.kernel, p0.dis_sampler.s
for x in (0.01, 0.10, 0.30):
    q2 = 5.0
    t = k.tables(x, q2)
    azz = lg.azz(t.b1, t.f1, t.f2, x, q2 / (x * s), t.b2, 0.0)   # 6Li b1
    rU  = m0.tail_ratio_at(x, q2, 0.0)
    rT0 = m0.tail_ratio_at(x, q2, 1.0) - rU
    rT1 = m1.tail_ratio_at(x, q2, 1.0) - m1.tail_ratio_at(x, q2, 0.0)
    print(x, rU, rT0, rT1,
          (azz + 2 * rT1) / (1 + rU) - (azz + 2 * rT0) / (1 + rU),
          m0.delta(x) * abs(azz))
```

The `phase_C_numbers.md` reproduction of §B3.2 is the same loop with
`lg.toy_b1(x, q2, t.f1, lg.B1Mode.Digitized)` in place of `t.b1`.

`build/lipolgen_tests -tc="B3*" -s` prints the same three "tensor tail ×"
factors on the **test** sampler grid (40 × 24, not the CLI's 101 × 77):
**×1.28237 / ×3.74748 / ×1077.5**. Every table above is on the production grid
(×1.2824 / ×3.7475 / ×1005.8); the x = 0.30 factor is the one that differs
materially, because there the ratio is `1 + σ^q_U/σ^el_U` ≈ 10³ and the
denominator is an interpolation of a tail eight decades below its neighbours.

## B3.6 What B3 leaves open

* **The polarised quasi-elastic tail is still not computed.** Nothing here is
  a calculation. A real answer needs a POLRAD Eq. (44) tensor analogue built
  on the deuteron-in-⁶Li spin structure — the alignment of the nucleon momentum
  distribution, not the coherent quadrupole — and it does not exist for A = 6.
* **Zhou *et al.* is still not in the tree.** The neglect rests on a paper this
  repository cannot read, and the sentence quoted from it via HERMES is about
  a deuteron.
* **The right scale is unknown.** §B3.1's "one to two orders of magnitude" is a
  proxy argument from a deuteron *elastic* tail, not a measurement of anything
  quasi-elastic.
* **`phase_C_numbers.md` §8.1c/§8.2's `A_zz`/`τ`/`w_hi`/`ΔA_zz` columns are
  wrong by ×3.253983** (§B3.2) and were **not** republished here.

---

# B4 — the tagged `rc_tail ≡ 1` is half a kinematic fact, and the other half is now written down

## B4.0 The question, and the answer

`rc_tail` is identically 1 on every tagged channel "because the tag itself
vetoes the elastic recoil at x_L = 1". **The elastic half of that is a fact.
The quasi-elastic half is false**, and this task worked out what is true
instead. Nothing was implemented; a **documented exclusion** replaced an
unexamined assumption, and the reason string, the header, four documents and
two test suites now carry it.

## B4.1 The kinematics, in order

1. **The elastic tail really is vetoed.** `e + A → e' + γ + A(g.s.)` leaves the
   ion **intact**, at `x_L = 1` and inside the 10σ beam-exclusion envelope
   (`PipelineConfig::n_sigma = 10`), while the tag looks at
   `x_L ≈ A_spec/A_beam` = 2/3 (⁶Li α), 4/7 (⁷Li α), 1/2 (d control). Not by
   rigidity — ⁶Li and α share `Z/A = 0.5` — but by **longitudinal fraction**.
   That is a property of the route classification and is not an approximation.
2. **A quasi-elastic knockout does not leave the ion intact.** It removes one
   nucleon. The A−1 remnant of either lithium channel is **unbound**: ⁵Li and
   ⁵He have no particle-stable state at all, both being resonances above the
   α + N threshold. So the remnant breaks up, and its α carries
   `x_L ≈ (4/5)(5/6) = 2/3` — **the tag window itself** — with a `p_T` spread
   of the order of the Fermi motion the tag already accepts.
3. **In the α + d picture this generator actually uses, it is sharper than
   that.** If the struck nucleon is one of the embedded **deuteron's** two, the
   α is a **true spectator**: nothing touches it, and its momentum distribution
   is the same `n_M(k, c)` the tagged Born is built on. On the deuteron control
   the statement is bare — the quasi-elastic tail there **is** elastic `e`–`n`
   scattering with a spectator proton, the classic spectator-tagging
   background.
4. **So the tag does not even suppress it in the ratio.** `rc_tail` is a ratio.
   The same spectator density and the same Roman-Pot acceptance `ε(k, c)`
   multiply its numerator and its tagged Born, and **both** select the 2-of-6
   nucleons inside the deuteron: a nucleon knocked out of the **α** destroys
   the α and *is* vetoed, which is the other 4-of-6. To leading order those
   factors cancel, and the omitted dilution is of the **same order as the
   inclusive quasi-elastic one** — 22 % / 73 % / 99.9 % of the inclusive tail
   at x = 0.01 / 0.10 / 0.30.

**That is an order-of-magnitude argument, not a computed number**, and it is
labelled as one everywhere it appears. What it refutes is specific: the
design's *"it is probably small (a radiative tail of a quasi-elastic peak into
a DIS bin, with the α-tag's own acceptance on top)"*. The α acceptance is not
"on top" — it is **common to numerator and denominator** and cancels.

## B4.2 Why nothing was implemented

Three things are missing, and none is a small piece of code:

* **A tagged Born denominator.** `RcModel::born_pb_at` is the **inclusive**
  `x·s·dσ_unpol`. The tagged Born is the struck-cluster DIS cross section times
  the spectator density and the pot acceptance, which lives in `TaggedModel`
  and is not reachable from the tail tables.
* **The tag acceptance inside Eq. (44).** POLRAD Eq. (44) is an inclusive
  incoherent sum over Z + N free nucleons. Restricting it to "the struck
  nucleon was in the deuteron **and** the α survived **and** landed in the
  pots" is a different integral, not a factor.
* **A piece Eq. (44) does not contain at all**: the **cluster-elastic**
  `e + A → e' + γ + d + α`, i.e. elastic scattering off the embedded deuteron
  as a cluster with the α as spectator. Its final state **is** the tagged one,
  its tag acceptance would be perfect, and its radiative tail reaches down to
  `Q'²` where `F_d` is alive. It is not in the free-nucleon sum and it is not
  priced anywhere.

Guessing any of these would have been a second uncited model layered on the
first. The exclusion stays, and the run says what is excluded.

## B4.3 What changed, and what gates it

`RcModel::exclusion_reason()` on every `Tagged*` channel now states the FACT
and the OMISSION separately, names the `x_L ≈ 2/3` α, the true-spectator
argument, the 2-of-6 / 4-of-6 split, the 22 / 73 / 99.9 % order of magnitude
and the cluster-elastic piece. It travels into the run banner **and** into
`meta["rc_exclusion_reason"]`, which is where an analysis reading
`rc_tail == 1` out of an npz would otherwise stop.

The word **`ASSUMPTION` was removed**, and a doctest asserts it is *absent*:
calling the quasi-elastic half an unquantified assumption is exactly what this
task replaced, and the word must not come back without the argument. The
positive tokens (`x_L`, `FACT`, `OMISSION`, `QUASI-elastic`, `2/3`,
`SPECTATOR`) are asserted present in both suites — `tests/test_rc.cpp` T4 and
`python/tests/test_rc.py::test_tagged_tail_exclusion_names_the_quasi_elastic_omission`,
which also checks the string reaches `meta`.

No weight moved: `rc_tail` is still exactly 1.0 on every tagged event, the band
is untouched, and the tagged clip statistics are unchanged.

## B4.4 What B4 leaves open

* **Q3 is still open**, and now with a *larger* estimate than it had. The
  quasi-elastic tagged tail is unpriced, and the argument above says it is not
  small.
* **The cluster-elastic tail is unpriced on every channel**, tagged or not —
  Eq. (44) is a free-nucleon sum and nothing in the tree computes coherent
  scattering off the embedded deuteron.
* **The "2 of 6" counting is a light-cone-fraction argument, not a
  spectroscopic one.** It ignores the α's own internal Pauli blocking, the
  difference between the struck-nucleon and struck-deuteron light-cone
  kinematics, and `S_αd`. Each is an O(1) factor; together they are why the
  claim is "same order", not a number.

## B4.5 Suite tallies

| | before B3/B4 (B1 + B2 applied) | after B3 + B4 |
|---|---|---|
| `build/lipolgen_tests` | 375 cases / 17 215 875 assertions / 1 skipped / 0 failed | **376** cases / **17 215 977** assertions / 1 skipped / **0 failed** |
| `python -m pytest python/tests -q` | 264 passed | **267 passed** |
| `validation/check_physics_channels_links.py` | 1061 refs / 0 broken | **1064 refs / 0 broken** |

One new doctest case (B3) and 102 new assertions; B4's are inside the existing
T4 case. Three new pytests (the B3 chain, the B3 refusals, the B4 exclusion
text) and three existing ones extended.

`--fix` was run (this task moved ~140 lines in `rc.hpp` and ~55 in `rc.cpp`)
and **its diff was reviewed mechanically, not trusted**. Of 130 references
whose line moved, the fixer re-pointed **11** from a declaration to a *comment
mention* of the same symbol — the same resolver failure §B2.7a documented — and
all 11 were repaired by hand — ten of them had landed on a comment, and
`delta` had landed on a *different declaration's* trailing comment
(`RcWeights::lo`, whose `///<` mentions `delta(x)`), which the comment scan
alone would have missed. Final targets after a second `--fix` pass (a later
header edit shifted `rc.hpp` again by five lines): `RC_X_HIGH` → `rc.hpp:192`,
`delta` → `:1228`, `LI6_R2_POINT_FM2` → `:503`, `applies` → `:1180`,
`clipped_cell_fraction` → `:1329`, `weights` → `:1333` and `rc.cpp:1621`,
`fill` → `:1340`, `RcModel` → `:1128`, `tail_ratio_at` → `rc.cpp:1362`,
`lam_e` → `rc.cpp:975`.

Verified **mechanically**, not by eye: the pre-edit `rc.hpp`, `rc.cpp` and
`pipeline.cpp` were reconstructed by reverse-applying this task's own text
swaps, and the **text** of every moved reference's target was compared before
and after. Of **97** moved references into those three files, **0** point at
different text than they did (the 98th is this task's own new
`qe_tensor_scale` reference, which resolves to the declaration). Eight
references — five distinct symbols: `make_tagged`, `default_inclusive_kernel`,
`beta` and `hadronizer` in `pipeline.cpp`, `LI6_VMC_C0_Q_MATCH_FM` in `rc.hpp`
— were **already** pointing at comment mentions before this task ran and were
left as found: they are not this task's regressions and repairing them was not
in scope.


---

# B5 — a transcription error that was only ever in the document

## B5.0 The finding, and the order it was established in

`polrad_transcription_check.md` §7.1 recorded that
`design_C_tensor_rc.md` §1.4.6 printed POLRAD Eq. (A.4)'s `ℑ^el_6` with a
**spurious `η_A`** in the `F_q` numerator. The brief's instruction was to
confirm the CODE is right and the DOCUMENT wrong **before** changing anything,
and to stop and report if it turned out the other way. It did not: the code is
right, and by the strongest possible route — **it does not contain `ℑ^el_6`
at all.**

**1. The source, read directly.** `polrad2t.tex` lines 2705–2706, verbatim:

```tex
\Im ^{el}_{6}= {Q_N\over 24} \pmatrix{F_m^{2}
    +{4\over 1+\eta_A}({\eta_A\over 3}F_q+F_c+\eta_A F_m)F_q} ,
```

`4/(1+η_A)`. No `η_A` in that numerator. (`ℑ^el_2`'s neighbouring
`4η_A²/(1+η_A)` **is** correct as the design prints it — the two are not the
same coefficient, which is presumably how the slip happened.)

**2. The FORTRAN, independently.** `ℑ₆ = (Q_N/6)ε³(b₂/3 + b₃ + b₄)` (Eq. (A.3))
with `strf`'s elastic `b₂, b₃, b₄` (adgh:4235–4242) and `ε = 1/(2τ)`. Carried
out term by term:

* `f_m²`: `(4/3)τ²[1 − (3τ+2) + (6τ+1)] = (4/3)τ²·3τ = 4τ³`
* `f_q`: `(16/9)(τ³/τ₁)[(τF_q+3F_c−3F_m) + (τF_q+3F_c−3F_m) + (τF_q+3F_c+3(3τ+2)F_m)]F_q`
  `= (16/9)(τ³/τ₁)[3τF_q + 9F_c + 9τF_m]F_q = (16/3)(τ³/τ₁)(τF_q+3F_c+3τF_m)F_q`

so `ε³(b₂/3+b₃+b₄) = (1/2)[F_m² + (4/(1+η))((η/3)F_q + F_c + ηF_m)F_q] · 2M_A²`
— **`4/(1+η_A)`**, with no second `η_A` anywhere, in POLRAD's own
`ℑ_code = 2M_A²·ℑ_paper` normalisation.

**3. The code.** `RcTailModel::PolradFull` — the only path on which Eq. (A.4)
at `Q_N ≠ 0` would ever be evaluated — throws
`"NOT implemented in v0"` in `RcModel`'s constructor
(`src/core/rc.cpp:933-935` (as of d3ac125)). What v0 actually evaluates is Eq. (38)
(`polrad_sigma_el_u` and the `σ_q` integrand, `src/core/rc.cpp:258-294`) and
the `Q_N = 0` Rosenbluth pair (`rosenbluth_spin1`, `src/core/rc.cpp:268-272`).
A `grep` for the coefficient finds nothing: **there is no `ℑ^el_6` in the
tree to be wrong.** So this is a **document-only** defect, and the larger
finding the brief warned about does not exist.

## B5.1 What was changed

Two files, both annotated rather than rewritten:

* `docs/open_items/run_2026-09-02/design_C_tensor_rc.md` §1.4.6 — the one line
  inside the Eq. (A.4) block now reads `(4/(1 + eta_A))`, and a dated
  **CORRECTED 2026-09-04 (task B5)** box directly under the block prints what
  the document said before, both independent verifications above, and the fact
  that the code never carried it.
* `docs/open_items/run_2026-09-02/polrad_transcription_check.md` §7.1 — an
  **APPLIED 2026-09-04** box, so the check that found it records that it was
  acted on and what was re-verified first.

**Nothing else moved.** No source file, no test, no number. The suites and the
docs gate were run anyway (§B6.6) because the same commit carries B6.

---

# B6 — the A = 2 → A = 6 transfer of the band gets a price, and the Gakh–Shekhovtsova SHAPE is refused on the record

## B6.0 What this task did, and what it deliberately did not

`δ(x)` is the **dominant** RC systematic — at `x = 0.01, Q² = 5` the band
half-width is `1.358472e−04` against the whole tail's
`ΔA_zz = −1.769797e−07`, a factor **768** (see §B6.5 on why the number quoted
elsewhere is 1300). Every anchor `rc_delta` interpolates between is a
**deuteron** number, and design question **Q8** calls whether the deuteron's
*fractional* RC transfers to ⁶Li "the largest unquantified assumption in the
band". Until now that assumption had a sentence and no dial: the shipped band
**asserted** `f = 0` — perfect transfer — silently.

**Priced, not corrected.** `RcOptions::a_transfer_frac` (default **0.0**) is
now that dial. **No A > 2 tensor RC calculation exists**, so none of `0`,
`0.5`, `1` is derived from one and the document says so in those words.

**And a rejection was recorded rather than left to be re-discovered.** The
obvious next move — replace the log-linear `rc_delta` with the shape from the
source's own figures — is refused, with its numbers, in `design_C_tensor_rc.md`
Q7 and `OPEN_ITEMS_SOLUTIONS.md` §9. The digitisation was **not** attempted;
what was done is the measurement that decides it (§B6.3).

## B6.1 The knob

```cpp
// include/lipolgen/rc.hpp, RcOptions, the "--- the band ---" block
double a_transfer_frac = 0.0;
```

read at exactly one place — `RcModel::delta` (`src/core/rc.cpp:1639`), which
is the **whole model's only call site of `rc_delta`**, so the term cannot be
applied twice or skipped on a path:

```cpp
const double d = rc_delta(x, opt_.delta_high_x, opt_.delta_low_x,
                          opt_.x_high, opt_.x_low);
return std::hypot(d, opt_.a_transfer_frac * d);
```

i.e. `δ_eff(x) = δ(x)·√(1 + f²)`. Three properties, each deliberate:

1. **Quadrature**, because it is an independent uncertainty on the same
   quantity, not a second correction to add.
2. **Multiplicative and δ-proportional**, so it vanishes where `δ` does. The
   E12-13-011 high-`x` anchor is an `A = 2` measurement too; the doubt about
   the transfer cannot be *larger* than the correction it doubts.
3. **`rc_delta` itself is untouched.** The free function stays the published
   shape — what T3 gates with `==` on both anchors and what the Python
   `rc_delta` binding exposes. Only the model's own `delta()` carries the
   price.

A negative `f` is **refused** (`RcModel`'s constructor and
`PipelineConfig::validate()`, both matched by tests) because the quadrature
erases its sign and `meta` would then record a price that was not charged.

## B6.2 The measurement the brief asked for

⁶Li, `--config 1` (`s = 3980 GeV²`), `--channel inclusive`,
`--plan tensor-thirds --pzz 0.6`, defaults, `ho`, `θ_S = 0`, **x = 0.01,
Q² = 5 GeV²**. `A_zz = −4.528242e−04` — this generator's own ⁶Li b₁, i.e. the
corrected reference of §B3.2 and **not** `phase_C_numbers.md` §8.2's
×3.253983-too-large deuteron-b₁ value. Half-width on `A_zz` = `δ_eff(x)·|A_zz|`:

| `a_transfer_frac` | `δ_eff(0.01)` | **band half-width on `A_zz`** | widening | band ÷ whole tail |
|---|---|---|---|---|
| **0 (default)** | 0.30000000 | **1.358472e−04** | 1× | ×768 |
| 0.5 | 0.33541020 | **1.518818e−04** | 1.118034× | ×858 |
| 1 | 0.42426407 | **1.921170e−04** | 1.414214× | ×1086 |

The widening is `√(1+f²)` to the last digit by construction. At the other two
published `x`, same `Q²`: `9.680016e−05 → 1.082259e−04 → 1.368961e−04` at
`x = 0.10`, and `7.307492e−07 → 8.170024e−07 → 1.033435e−06` at `x = 0.30`.

**Reading.** `f = 1` does not change the ordering of the RC budget — the band
still leads the whole tail by ×1086 and every tail knob by more — it changes
the *size* of the thing that already dominates, by 41 %. That is the honest
shape of this systematic: the largest term in the budget is also the one with
no `A = 6` support at all.

## B6.3 Why the Gakh–Shekhovtsova SHAPE is NOT adopted — measured, and recorded

The four `(x, Q², δ)` panels of hep-ph/0403262 Fig. 2 are available as **data**,
not as a picture: the arXiv source (`https://arxiv.org/e-print/hep-ph/0403262`)
ships `rrc01.eps`, `rrc1.eps`, `rrc4.eps`, `rrc10.eps` for
`Q² = 0.1, 1, 4, 10 GeV²`, each with a 30-point `/x` array and **two** `/y`
arrays — the first drawn under `gsave [25 10] 0 setdash` (the dotted,
RC-included curve of the caption), the second solid (Born). Everything below is
read off those arrays.

**The `δ` used below, and why it is that one — CORRECTED 2026-09-04.**
`δ ≡ (Δσ_RC − Δσ_Born)/Δσ_Born`, the radiative correction as a fraction **of
the Born**. The paper never writes a `δ` at all — it plots two curves and quotes
one percentage range — so the *ratio* is this record's own construction:
**inferred, not sourced**. What **is** sourced is the DENOMINATOR, and it decides
the form. The paper's one relevant quoted number is *"the value of radiative
correction changes from 10 % to 30 % **as compared with the Born contribution**"*
(`DIS-TENS2.tex` line 924), which names `Δσ_Born` as the base — and that is the
base under which it reproduces. Measured on panel (a) — `Q² = 0.1`,
`x = 0.00226 – 0.00966`, the **only** panel lying inside the paper's own quoted
`x ∼ 10⁻³–10⁻²`:

| definition of `δ` | `\|δ\|` across panel (a) | to 1 significant figure |
|---|---|---|
| **`(Δσ_RC − Δσ_Born)/Δσ_Born`** — ADOPTED | **0.1133 … 0.2662** | **10 % … 30 %** — the paper's sentence |
| `Δσ_Born/Δσ_RC − 1` — this section's original | 0.1278 … 0.3628 | 10 % … **40 %** |

**The earlier version of this section asserted the opposite** — that
`Δσ_Born/Δσ_RC − 1` was "the definition under which the paper's own quoted
numbers reproduce". **That attribution was false.** It was inferred by matching
the task brief's values, never read from the paper, and the paper's number
reproduces *better* under the form the section did not adopt. The `δ_low = 0.30`
anchor is independent confirmation: it was taken from the paper's "30 %", and it
is the adopted `δ` whose panel maximum — 0.2662 — rounds there; the original
definition's 0.3628 would have set that anchor at 0.40.

Two further reasons the adopted form is the right one for **this** repository,
independent of the paper: (i) the band ansatz is `w = 1 ∓ δ·τ` applied to the
generator's **Born** rate, so `δ` must be a fraction *of the Born* — the original
form puts the *corrected* rate in the denominator, which is not the quantity the
weight multiplies; (ii) `0.113…0.266` was already published under the symbol
`δ` in the shipped `OPEN_ITEMS_SOLUTIONS.md` §9 honest flags, so adopting it
leaves **one** definition across every file, as `docs/CONVENTIONS.md` requires.
Everything below is quoted under the adopted `δ` **only**; a negative `δ` means
the RC *reduces* `|Δσ|`.

**(a) δ is SINGULAR, not a bounded fraction.** The Born `Δσ` crosses zero — the
`b₁` zero-crossing — inside the `Q² = 4` and `Q² = 10` panels, and the RC
**moves the crossing**, which is the paper's own stated result (*"the inclusion
of the radiative correction shifts the zero value of `b₁` and `b₂` to the
smaller `x`-value region (see Fig. 2 c and d)"*). Measured:

| `Q²` | `x` window | Born `x₀` | RC `x₀` | shift | **δ range on the panel** |
|---|---|---|---|---|---|
| 0.1 | 0.00226 – 0.00966 | — | — | — | −0.266 … −0.113 |
| 1 | 0.01348 – 0.09662 | — | — | — | **−0.171 … +0.240** (sign change) |
| 4 | 0.05394 – 0.38647 | 0.20118 | 0.18467 | −0.01652 | **−1.691 … +7.810** |
| 10 | 0.13097 – 0.85000 | 0.20161 | 0.18331 | −0.01830 | **−0.887 … +4.091** |

The `x₀` columns are properties of the two curves and do not depend on the
definition; the `δ` column does, and is the adopted one. Under it the **pole sits
at the Born zero** (the denominator's), and the panel extrema are the two grid
points that straddle it: −1.691 at `x = 0.19154` and +7.810 at `x = 0.20300`
(`Q² = 4`, spacing 0.01147, `x₀ = 0.20118`); −0.887 at `x = 0.18056` and +4.091
at `x = 0.20535` (`Q² = 10`, spacing 0.02479, `x₀ = 0.20161`). On panel (b) the
sign change runs **+0.240 at the low-`x` end to −0.171 at the high**.

At `Q² = 4` in detail: **−0.170** (`x = 0.111`) → **−0.537** (`x = 0.169`) →
**+1.026** (`x = 0.214`) → **+0.030** (`x = 0.283`). The band ansatz is
`w = 1 ± δτ` with `band_tau_max = 1`, which assumes a small, bounded,
one-signed fractional rescale and needs `|δ| < 1` merely to keep both weights
positive. A
paper-shaped `δ` is unbounded, changes sign, and would emit **negative
weights**.

**(b) It is not monotone non-increasing in `x`.** That monotonicity is the
design's stated justification for the log-linear form (§1.3) and is *asserted*
by `tests/test_rc.cpp:70-83` — T3 walks a 2000-point log grid and compares
**both anchors with `==`**. Beyond the sign changes of (a), the panels give
**−0.128 at `x = 0.386`** (`Q² = 4`) and **−0.400 at `x = 0.85`** (`Q² = 10`)
against `δ_high = 0.015` for **all** `x ≥ 0.16` — **8.5×** and **26.7×** past the
E12-13-011 anchor in magnitude, and negative besides, which a half-width
cannot be. The two anchors and the shape cannot both stand; adopting
the shape would silently delete the only *measured* statement in the band.

**(c) Patchy, non-rectangular `Q²` support.** Four `Q²` values whose `x`
windows barely overlap, and panels **a and b do not overlap at all** — there is
a gap at `0.00966 < x < 0.01348`. This generator's own quoted point,
**`x = 0.01, Q² = 5`**, falls *between* the `Q² = 4` and `Q² = 10` panels in
`Q²` **and** outside the `x` window of both (`0.01 < 0.05394`,
`0.01 < 0.13097`), with its `x` landing exactly in the a–b gap. **There is no
panel to interpolate from at the point the design publishes.**

**(d) Their Born is not our Born.** The ratio is taken against *their* tensor
model: Eq. (60)'s HERMES parametrisation
`A_zz = −1.56·10⁻²(1 − 1.74x − 1.45√x)`, `F₂^d` from ALLM97 [21] with
`F₂ⁿ/F₂ᵖ` [22], and `b₂` from the Callan–Gross-type relation Eq. (59) — all
four read verbatim in the paper's own §4 *Numerical estimations*, not inferred.
(This section said "§5" until 2026-09-04; the arXiv source has four numbered
sections plus three appendices, and *Numerical estimations* is the fourth. The
attribution was right, its section number wrong. It is **not** a precedent for
the `δ` definition above, which is inferred and now says so.) Transferring a
*fractional* RC computed against a different `b₁` shape onto this generator's
`b₁` is a fresh unpriced assumption stacked on top of Q8's.

**A by-product worth its own flag — CORRECTED 2026-09-04, twice over.** The
*magnitude* the design does take, the paper's "10 % to 30 %", is stated for
`x ∼ 10⁻³–10⁻²`, and **exactly one panel lies in that window**: panel (a),
`Q² = 0.1`, `x = 0.00226 – 0.00966`, all 30 of its points inside it. Panel (b)
starts at `x = 0.01348`, *above* `10⁻²`, and has **no** point in the quoted range.
So `δ_low = 0.30` is a **`Q² = 0.1 GeV²`** number quoted at `Q² = 5` — one panel,
not the "`Q² = 0.1–1 GeV²`" two this section claimed until now, which makes the
flag *stronger*. On panel (a), `|δ| = 0.113 … 0.266`.

**And that pair is an `x` spread, not a `Q²` spread — with a consequence this
record did not state until now.** Panel (a) being the only panel in range, the
"10 %" and the "30 %" are its two **`x` ends at one `Q²`**: `|δ| = 0.1133` at
`x = 0.00966`, the panel's top, and `0.2662` at `x = 0.00226`, its bottom. The
paper's own sentence scopes them the same way — *"In the range of **low x**
(x ∼ 10⁻³–10⁻²) … changes from 10 % to 30 %"* — and names no `Q²` at all. Four
shipped sites nonetheless described the pair as *"the spread of ONE calculation
over Q², not two independent edges"*: `PHYSICS_CHANNELS.md`'s `δ(x)` row and its
`[GS04]` entry, `include/lipolgen/rc.hpp` immediately above `RC_DELTA_LOW_X`,
and `design_C_tensor_rc.md`'s low-`x` bullet. **All four are corrected
2026-09-04.**

**THE CONSEQUENCE, which no file stated before this one.** `δ_low = 0.30` is
anchored at `RC_X_LOW = 0.01`, and 0.01 lies **above** panel (a)'s top
`x = 0.00966`, where the panel itself reads `|δ| = 0.113`. The `0.266` that
rounds to the paper's "30 %" is the panel's value at `x = 0.00226` — a factor
**4.3 lower in `x`**. So the anchor's real basis is **"the panel's lowest-`x`
value, carried upward in `x` to 0.01"**, which is *not* the same statement as
**"the conservative end of a `Q²` spread"** — and the shipped documents
asserted the second. Both readings are conservative in *magnitude*
(`0.30 > 0.266 > 0.113`), and magnitude is the only property `w = 1 ∓ δτ` uses,
so **nothing moves**: this is recorded, not fixed. Re-siting an anchor after
seeing why it sat where it did is the author's call, not a correction pass's.

**And `0.171…0.240` was not a `|δ|` range under either definition.** Under the
adopted `δ`, `−0.171` and `+0.240` are panel (b)'s two **signed extremes**, at
opposite ends of its `x` window; `|δ|` on that panel actually spans
**0.006 … 0.240**, so the published interval excluded most of the panel and its
lower "edge" was a *sign*, not a magnitude. Both errors are corrected here and
in the **two** other files that carried them: `OPEN_ITEMS_SOLUTIONS.md` §9 and
`design_C_tensor_rc.md` Q7. (**Corrected 2026-09-04**: this sentence used to say
"the three other files that carried them" and name `STATUS.md` row 7 as the
third. `STATUS.md`'s decision row 7 is not a third — it was written in this same
phase, after the correction, and never carried either the `0.171…0.240` range or
the `Q² = 0.1–1 GeV²` claim; it states the one-panel reading directly.)
Recorded in
`OPEN_ITEMS_SOLUTIONS.md` §9's honest flags; it makes the case against the
*shape* stronger, not weaker, since taking the magnitude is already one
extrapolation.

**What is kept:** the log-linear interpolation, the two visible anchors, and the
magnitude with its flags. Q7 stays **open**; the digitisation route is now
**closed**, and reopening it needs a `δ` definition that survives a zero
crossing.

## B6.4 The `δ_low` source still has ZERO INSPIRE citations — re-verified

Checked 2026-09-04 against the INSPIRE literature API
(`inspirehep.net/api/literature?q=arxiv:hep-ph/0403262`): record **647050**,
*"Radiative corrections to deep inelastic e d scattering: Case of tensor
polarized deuteron"*, **`citation_count: 0`**. **The flag stands** and is kept
wherever it was, in `rc.hpp`'s `RC_DELTA_LOW_X` comment, the run banner,
`USAGE.md`, `PHYSICS_CHANNELS.md` and `OPEN_ITEMS_SOLUTIONS.md` §9.

## B6.5 The default did not move — measured, not argued

`std::hypot(d, 0.0)` returns `|d|` exactly and `rc_delta` is non-negative at
any sane anchors, so the shipped band is unchanged *in principle*. It was
checked *in practice*, by the §B3.4 method: the one-line change was reverted in
place, the library rebuilt, two CLI runs regenerated at the same seed, and the
files compared column by column as `ndarray.tobytes()` plus a key-by-key `meta`
diff.

```
python -m lipolgen.cli --isotope 6Li --channel inclusive --config 1 \
    --events 2000 --seed 4242 --rc {tensor-band,off} --npz …
```

| run | columns compared | byte-differing | `meta` keys | differing |
|---|---|---|---|---|
| `--rc tensor-band` | **51** | **0** | 51 | **0** |
| `--rc off` | **48** | **0** | 24 | **0** |

`--rc off` was additionally compared against a run made **before** this task
(the B3-era build): 48 columns byte-identical, 24 `meta` keys identical, none
added or removed. The `--rc off` guarantee is intact.

**And the knob does exactly one thing when turned on.** `f = 0` vs `f = 1`,
same seed:

| | value |
|---|---|
| columns that differ | **`rc_tensor_lo`, `rc_tensor_hi` — and no others** |
| columns byte-identical | **49** (including `rc_tail`, `weight`, every kinematic) |
| `meta` keys that differ | **`rc_a_transfer_frac`** only (0.0 → 1.0) |
| `(w_hi − 1)` ratio, mean over 2000 events | **1.4142135623730436**, max deviation from `√2` **5.3e−11** |

The `5.3e−11` is the cancellation in `w_hi − 1`, not a modelling error:
`w_hi` sits within `1e−4` of 1, so the subtraction carries `~1e−16/|w_hi−1|`
of relative noise. The same effect sets the C++ test's `1e-9` tolerance on that
identity, and the comment there says so.

**On the factor 1300.** The brief quotes the band as "a factor 1300 larger than
the whole tail at `x = 0.01, Q² = 5`". That number, like the `4.4e−04`
half-width it came from, was built on the deuteron-b₁ `A_zz` corrected in
§B3.2. With this generator's ⁶Li b₁ the ratio is **768** (`1.358472e−04` ÷
`1.769797e−07`). **The ordering the row exists to state is unchanged** — the
band still dominates the tail by nearly three orders of magnitude — and the
arithmetic is now right at both sites (`phase_C_numbers.md` §8.1 item 3 and
`OPEN_ITEMS_SOLUTIONS.md` §9 item 3, each annotated in place).

## B6.6 What gates it, and how to reproduce

**C++**, `tests/test_rc.cpp`, `TEST_CASE("B6: the A = 2 -> A = 6 transfer
price…")`, four subcases: the default is 0 and `RcModel::delta` **`==`**
`rc_delta` on a ten-point `x` grid including both anchors; `δ` widens by
exactly `√(1+f²)` at `f = 0.25, 0.5, 1, 2` while the free `rc_delta` is
unmoved; the knob leaves `tail_ratio_at` and `tensor_fraction` **bit-identical**
across the sampler grid while scaling `w_hi − 1` and `w_lo − 1` by `√2`; and a
negative `f` throws.

**Python**, `python/tests/test_rc.py`,
`test_a_transfer_frac_prices_the_A2_to_A6_transfer_of_the_band` — the `meta`
key exists and is 0.0 by default, the round trip through `make_config` and the
CLI default works, `delta` widens by `√(1+f²)` at `f = 0.5, 1, 2` with the
tail untouched, and a negative `f` is refused by `validate()`. The `meta`
key-set assertion in `test_meta_records_every_rc_knob_that_moves_a_column` now includes
`rc_a_transfer_frac`, which is what stops the knob from ever becoming
invisible in an npz.

```bash
source env.sh
cmake --build build -j
build/lipolgen_tests -tc="B6*"          # 1 case, 649 assertions
python -m pytest python/tests/test_rc.py -q
python3 validation/check_physics_channels_links.py
```

The §B6.2 table, from the shipped Python API alone:

```python
import lipolgen as lg
from lipolgen import _lipolgen as _l
def mk(**kw):
    cfg = lg.make_config(isotope="6Li", channel="inclusive", config=1,
                         events=2000, seed=99, rc="tensor-band", **kw)
    return lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
p0 = mk(); k, s = p0.dis_sampler.kernel, p0.dis_sampler.s
x, q2 = 0.01, 5.0
t = k.tables(x, q2)
azz = lg.azz(t.b1, t.f1, t.f2, x, q2 / (x * s), t.b2, 0.0)   # 6Li b1
for f in (0.0, 0.5, 1.0):
    m = (mk(rc_a_transfer_frac=f) if f else p0).rc_model
    print(f, m.delta(x), m.delta(x) * abs(azz))
```

The §B6.3 panel table, from the arXiv source of the paper:

```bash
curl -sL -o src.tar.gz https://arxiv.org/e-print/hep-ph/0403262
tar xzf src.tar.gz            # rrc01.eps rrc1.eps rrc4.eps rrc10.eps
```

then, per file, parse the `/x [...] def` array and the **two** `/y [...] def`
arrays (first = dotted = RC, second = solid = Born) and form
`(rc[i] − born[i])/born[i]` — the adopted `δ` of §B6.3. (Until 2026-09-04 this
recipe said `born[i]/rc[i] − 1`; it now matches the tables it reproduces.)

## B6.7 Suite tallies

| gate | before Phase B (`ac22331`) | after B5 + B6 |
|---|---|---|
| `build/lipolgen_tests` | 374 cases / 17 208 032 assertions / 1 skipped / 0 failed | **377 cases / 17 216 856 assertions / 1 skipped / 0 failed** |
| `python -m pytest python/tests -q` | 261 passed | **268 passed** |
| `validation/check_physics_channels_links.py` | 1049 refs / 0 broken | **1071 refs / 0 broken** |

**The assertion total was re-measured 2026-09-04 and corrected here, +230.**
This table and `STATUS.md` row 11 both read **17 216 626** until then — the
count taken *before* the B2 acceptance correction added subcase **T8(d)(i)** to
`tests/test_rc.cpp`, whose assertions nothing re-tallied afterwards. The case
count (377) and the other two gates were already right. Attributed by
**measurement, not arithmetic**: `build/lipolgen_tests -sce="*LOW-y corner*"`
reproduces the old **17 216 626** exactly, so T8(d)(i) is the whole of the 230.
The full suite was run three times on the corrected tree, identical every time,
and the only files this task changed are Markdown.

**What B6 itself contributes:** `build/lipolgen_tests -tce="B6*"` gives
**376 cases / 17 216 207 assertions** (re-measured 2026-09-04; it read
17 215 977 before T8(d)(i) existed), so B6's one new case still carries
**649** assertions and B1–B4 account for the other two cases above. B5 adds none — it
changed no code. The docs gate was re-run with `--fix` after the source edits
shifted line numbers, and the resulting diff was reviewed **mechanically**: 39
lines of `PHYSICS_CHANNELS.md` changed, **101** individual `path:line`
references renumbered, and **0** lines whose non-line-number content changed
(compared under a normalisation that blanks the numeric part of every
reference).

One of the 101, `src/core/rc.cpp` → `weights`, resolves to the
`// ---- the weights` **section comment** rather than to
`RcWeights RcModel::weights(...)` at `:1643`. It did so before this task too
(it was `:1621`, the same comment, and moved by exactly the three lines this
task added); it is one of the pre-existing comment-anchored references §B2.7
recorded and **left as found**, not a regression introduced here.

## B6.8 What B5 + B6 leave open

* **Q8 is not closed.** There is still no `A > 2` tensor RC calculation, and
  `a_transfer_frac` is a price tag with no measurement behind any of its
  values. The default remains the *assumption* `f = 0`; it is now stated in the
  run banner instead of being silent.
* **Q7 is not closed either** — only one route to closing it is. The shape of
  `δ(x)` is still a choice; what B6 establishes is that hep-ph/0403262 cannot
  supply a replacement for it.
* **`δ_low = 0.30` is a `Q² = 0.1 GeV²` number** quoted at `Q² = 5` — one
  panel, corrected 2026-09-04 from "`0.1–1`" (§B6.3). It is flagged, not fixed; fixing it needs a calculation with `Q²`
  support where this generator runs.
* **`phase_C_numbers.md` §8.1c/§8.2's `A_zz`/`τ`/`w_hi`/`ΔA_zz` columns are
  still not republished.** B6 corrected the two *derived* statements that
  depended on them (the band half-width and the factor 1300) at their sites;
  the tables themselves remain as §B3.2 found them, with their warning.
* **`ℑ^el_6` remains unexercised.** The B5 correction makes the design record
  right; nothing evaluates Eq. (A.4) at `Q_N ≠ 0` until
  `RcTailModel::PolradFull` is built, and T9's `PolradFull` half is still
  `SKIP`-ped.
