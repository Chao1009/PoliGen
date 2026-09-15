# Phase B numbers — 2026-09-06 follow-on run

## B0. Baseline re-measured, read-only

Measured fresh in this session, tree clean throughout, nothing built (all
targets already up to date — `cmake --build build -j` printed only "Built
target ..." lines, no compile/link step ran).

| check | command | result |
|---|---|---|
| C++ build | `cmake --build build -j` | no-op (all targets already built) |
| C++ tests | `build/lipolgen_tests` | **404 test cases / 17 240 391 assertions / 0 failed / 1 skipped** |
| pytest | `python -m pytest python/tests -q` | **1004 passed / 127 skipped** (161.51 s) |
| docs gate | `python3 validation/check_physics_channels_links.py` | main: **1220 references checked (strict), 96 ranges, 7 external, 0 broken, 6 allow-listed**; plus per-file: `SPIN32_FINITE_GAMMA.md` 19 checked / 115 ranges / 0 broken; `PYTHIA_BRIDGE.md` 0 checked / 8 external / 0 broken |
| SPDX gate | `python3 validation/check_spdx_headers.py` | **102 files checked, all carry `SPDX-License-Identifier: GPL-3.0-or-later`** |
| `git status --porcelain` | | empty |
| `git rev-parse --short HEAD` | | `53f6951` |

**Reconciling this against `STATUS.md`.** `STATUS.md`'s own baseline line (top
of the file) quotes the *2026-09-06 run's start* state, commit `5af0427`
(before Phase A and the CW fix). The numbers above are the *current* HEAD,
`53f6951` — one commit past the CW-fix status-board update (`1d7e92a`), adding
`53f6951` itself ("sign fix, residues: ... tightened 2e-3 -> 1e-4"), which
touches only five doc files (`docs/OPEN_ITEMS_SOLUTIONS.md`,
`docs/USAGE.md`, `docs/benchmarking/07_cw_sign_investigation.md`,
`docs/open_items/run_2026-09-06/phase_CW_numbers.md`,
`docs/open_items/vmc_reconciliation.md`) and does not touch `STATUS.md`, so
`STATUS.md`'s own CW-row prose (the "Amended 2026-09-06 (verification pass)"
paragraph) is now one commit stale relative to HEAD — it still reads
"+3.1e-9" / "1e-3 -> 2e-3" where `53f6951`'s commit message records the
corrected "+3.1e-6" / "2e-3 -> 1e-4". This note does not fix that (read-only
phase); it is flagged here so the next phase does not quote `STATUS.md`'s CW
paragraph as the final word without checking `53f6951`'s diff.

The doctest/pytest/gate/SPDX tallies above match the CW row's own
"End state" line in `STATUS.md` (404 / 17 240 391 / 1 skipped; 1004 pytest /
127 skipped; gate 1220/96/7/0/6; SPDX 102/102) digit for digit — i.e. `53f6951`
moved no test, no gate count and no SPDX count, consistent with it being a
docs-only correction commit. **Nothing regressed. The batch reproduces from
committed state.**

---

## B1. The tensor fraction of the s-/p-peaks — bounded, not computed, and the bound is EMPTY where it matters

`SUMMARY.md` item 5: *"the **tensor** fraction of the s-/p-peaks is unknown,
not zero"*. Phase B of the 2026-09-03 run promoted the leading-log s-/p-peaks
into `src/core/rc.cpp` (`ll_peaks_spin1`, `ll_peaks_qe`, columns `u_sp`/`qe_sp`
of `TailTriple`) and deliberately kept them **out** of the tensor numerator of
`RcModel::tail_ratio_at`, because POLRAD's **Eq. (38)** supplies no tensor s/p
peak and deriving one from the leading log here would be an uncited second
definition. That left the tensor fraction of that piece at **exactly zero** — a
choice, not a measurement.

*(Say which POLRAD. §B2 below establishes that the **unqualified** sentence
"POLRAD supplies no tensor s/p peak" is true of Eq. (38) and **false of the
paper**: Eq. (18) + Eq. (A.4) carries it and `RcTailModel::PolradFull` computes
it. That is why `sp_tensor_scale` is refused on `PolradFull` for the
**opposite** reason it is refused on `TPeak` — not because the term did not
run, but because it did.)*

It is now **bounded, not computed**, by `RcOptions::sp_tensor_scale` /
`--rc-sp-tensor-scale` (default `0.0`). **And the number that goes with the
word is the point of this section: at the three standard points the bound is
not small, it is EMPTY — bit-identical to zero — and the piece that actually
carries the collapse is bounded by nothing at all.** Everything below was
measured in this session on ⁶Li `--config 1`, `s = 3980.0000000000005 GeV²`,
production grid (101 × 77 nodes, 3051 accepted cells), `ho` C0 edge, defaults
otherwise.

### B1.1 What was built

One term in `tail_ratio_at`'s numerator, beside the `u_sp` column it scales:

```
  num = su + (q_n/6) ratio_t su + su_sp
      + (q_n/6) * sp_tensor_scale * ratio_t * su_sp            <-- NEW
      + qe_suppression * ( sqe + sqe_sp + (q_n/6) qe_tensor_scale ratio_t sqe )
```

`sp_tensor_scale = 1` is therefore exactly the sentence **"the s/p tensor
fraction equals the elastic t-peak's"**. It sits **outside** `qe_suppression`,
because `u_sp` is the **elastic** s+p column and that knob is documented as a
flat multiplier on the *quasi-elastic* tail — the mirror image of
`qe_tensor_scale`, which rides *inside* it (2026-09-03 §B3.3 departure 1).

**Reach rule.** `PipelineConfig::validate()` and `RcModel`'s constructor both
**REFUSE** a non-zero value when `tail_model != TPeakPlusLL`, following the
`with_qe_tail` precedent and the criterion stated once at `KnobProvenance`: a
scale on a term the run computes as identically zero names *a variation of a
piece that did not run*, and no label makes a priced systematic un-priced.
Under `TPeak` the `u_sp` table is identically zero, so the setting would be
bit-identical to the default while `meta` recorded a price. The
`knob_provenance` row keys on exactly that:

| run | `rc_sp_tensor_scale` row |
|---|---|
| `--rc off` | `not-read` — "not read at rc = off: no RcModel is built" |
| `--rc tensor-band` (shipped `t-peak`) | **`refused`** — "tail_model = t-peak computes no leading-log s-/p-peaks at all (the u_sp table is identically zero)…" |
| `--rc tensor-band --rc-tail-model t-peak+ll` | `read` — "LEADING-LOG s-/p-PEAK knob… a BOUND WITH NO DERIVATION" |
| any tagged / coherent channel | `refused` by the tail rule, as every other tail knob is |

Both cells are in `python/tests/test_knob_provenance.py`'s matrix (`0 -> 1 at
t-peak+ll` and `0 -> 1 at t-peak (the shipped tail)`); the first moves the
60-event output hash, so the `read` row is earned and not asserted.

### B1.2 The default did not move — measured, not argued

The numerator was **restructured**, so bit-identity was checked empirically
rather than reasoned about: all six edited source files (`rc.hpp`, `rc.cpp`,
`pipeline.cpp`, `bindings.cpp`, `__init__.py`, `cli.py`) were reverted in
place, the library rebuilt, three npz regenerated at the same seed, and
compared column-by-column as `ndarray.tobytes()` plus a key-by-key `meta`
diff. The post-revert restore was checked by md5 on all six.

```
python -m lipolgen.cli --isotope 6Li --channel inclusive --config 1 \
    --events 2000 --seed 4242 --rc {off,tensor-band} [--rc-tail-model t-peak+ll] --npz …
```

| run | numeric columns | `meta` |
|---|---|---|
| `--rc off` | **48 / 48 byte-identical** | +`rc_sp_tensor_scale` **absent** (no RcModel); `knob_provenance` 68 → 69 rows, **no existing row changed** |
| `--rc tensor-band` (`t-peak`) | **51 / 51 byte-identical** | +`rc_sp_tensor_scale = 0`; 68 → 69 rows, none changed |
| `--rc tensor-band --rc-tail-model t-peak+ll` | **51 / 51 byte-identical** | +`rc_sp_tensor_scale = 0`; 68 → 69 rows, none changed |

And the whole-run 200 000-event `--rc tensor-band` run at seed 1234 is
**bit-identical between the two builds** (`rc_tail` array compared as raw
bytes), mean `1.021778527` on both.

> **CORRECTION, 2026-09-15 — investigated, and it is the FILL PLAN.** This
> paragraph carried *"That mean is **not** the `1.021836305` published by the
> 2026-09-03 run; the difference is **pre-existing** … and is not investigated
> here."* It is not a build difference at all. This comparison ran the CLI at
> its **defaults, P_z = 0.7**; the 2026-09-03 run published its numbers at
> **`--pz 0`** (`run_2026-09-03/phase_B_numbers.md:836-839` tabulates both).
> Re-measured on the current build, 200 000 events at seed 1234, `t-peak`:
> **1.021778527** at `tensor_thirds_plan(0.7, 0.6)` and **1.021836305** at
> `tensor_thirds_plan(0.0, 0.6)` — the 2026-09-03 digits, exactly. The
> bit-identity claim above (same plan, two builds) is untouched; see §B2.5's
> correction box.

`git diff --stat validation/reference/` is **empty**: no rtol-1e-12 reference
gate moved.

### B1.3 Why `TPeakPlusLL` lowers the tensor fraction, measured

`t-peak+ll` adds the s+p to the **unpolarised** numerator alone, so the tensor
fraction `r_T/r_U` falls purely because the denominator grew. At Q² = 5 GeV²:

| x | `r_U` (t-peak) | `r_U` (t-peak+ll) | denominator × | `r_T/r_U` (t-peak) | `r_T/r_U` (t-peak+ll) | **collapse ×** |
|---|---|---|---|---|---|---|
| 0.01 | 5.182218e−04 | 7.835318e−04 | ×1.51196 | −3.972572e−04 | −2.627428e−04 | **×0.66139** |
| 0.10 | 1.319634e−06 | 4.146075e−04 | **×314.184** | +4.161551e−05 | +1.324560e−07 | **×0.0031829** |
| 0.30 | 6.861763e−08 | 1.038467e−03 | **×15134.1** | +1.302390e−05 | +8.605653e−10 | **×6.6076e−05** |

The tensor NUMERATOR is literally the same expression on both edges —
`(q_n/6)·ratio_t·su`, untouched by `tail_model` — and it measures the same:
`r_T` = −2.058674e−07 on both at x = 0.01, and 8.936690e−13 vs 8.936688e−13 at
x = 0.30, where the 7th-digit difference is the round-off of differencing two
tail ratios whose common size grew by ×15 134. So the collapse is **entirely**
the denominator. That is the omission this knob is meant to price: a factor
**314** and **15 134** of new, entirely tensor-blind denominator at x = 0.10
and 0.30.

### B1.4 THE PRICE — and where the bound is empty

`ΔA_zz = (A_zz + 2 r_T)/(1 + r_U) − A_zz`, so `Δ(ΔA_zz) = 2[r_T(s) − r_T(0)]/(1 + r_U)`
is independent of `A_zz` (the knob carries `Q_N`; `r_U` does not move at all —
checked exactly at every accepted cell). Band half-width is
`δ(x)·|A_zz|` on this tree's own ⁶Li b₁ (the 2026-09-03 §B3.2 correction).
All rows on `--rc-tail-model t-peak+ll`, Q² = 5:

| x | band half-width | **Δ(ΔA_zz) at sp = 0.5** | **at sp = 1** | **% of band** | `qe_tensor_scale = 1` for contrast | % of band |
|---|---|---|---|---|---|---|
| 0.01 | 1.358472e−04 | **0** | **0** | **0.000 %** | −1.161702e−07 | 0.0855 % |
| 0.10 | 9.680016e−05 | **0** | **0** | **0.000 %** | +3.016431e−10 | 0.00031 % |
| 0.30 | 7.307492e−07 | **0** | **0** | **0.000 %** | +1.794093e−09 | 0.246 % |

**The zeros are exact, and they are the finding.** `sp_tensor_scale = 1`
returns a tail ratio **equal bit for bit** to `sp_tensor_scale = 0` at all
three points, for `q_n` = 0, +1 and −2. The reason is `RcOptions::sp_tensor_scale`
point (3), and it is measured rather than argued:

| x (Q² = 5) | y | t-peak vertex `t_min` | `F_c` there | s-peak vertex `Q'²_s = z_s Q²` | `F_c` there | `u_sp/σ^el_U` |
|---|---|---|---|---|---|---|
| 0.01 | 0.125628 | 8.7304e−05 GeV² | **+2.992599** | **4.37277 GeV²** | **−3.753238e−45** | 5.2492e−79 |
| 0.10 | 0.0125628 | — | — | **4.93822 GeV²** | **−6.435773e−51** | 2.0652e−86 |
| 0.30 | 0.0041876 | — | — | **4.98010 GeV²** | **−2.407836e−51** | 4.0603e−83 |

*(the `u_sp/σ^el_U` column is on the doctest's own smaller test-sampler grid,
the only place `TailTriple` is reachable — `tail_sigma_at` is not bound into
Python; the `F_c` and `Q'²` columns are grid-independent.
The p-peak is worse still: `Q'²_p = Q²/z_p = 5.717 GeV²` at x = 0.01, and at
the alive cell below it is 52.16 GeV², where `F_c` returns exactly −0.0.)*

The t-peak's ratio `σ^el_T/σ^el_U` is a quantity **measured at t ~ t_min**,
where the coherent form factor is alive. The s-/p-peak's own elastic vertex
sits **four decades higher**, at the scale of Q² itself, where it is **45 to 51
decades down**. So the borrowed ratio is not merely the wrong reference at high
Q² — the term it multiplies has already vanished, and the added tensor piece
falls below the numerator's ulp. **The bound is not conservative at the
standard points; it is empty there.**

**Where it is NOT empty.** The term is non-zero on **77 of 3051 accepted
cells**, every one with `x ≤ 7.943e−03` and `y ≥ 0.3656` (Q² from 1.054 to
30.91) — the low-x, `y → 1` corner, where `z_s = (1−y)/(1−x_A y) → 0` drags the
elastic vertex back down to where the form factor lives. There it is large:

| cell | `Q'²_s` | `F_c(Q'²_s)` | band half-width | **Δ(ΔA_zz), sp = 0.5** | **sp = 1** | **% of band** | `qe = 1` there |
|---|---|---|---|---|---|---|---|
| x = 4.16869e−04, Q² = 1.6081, y = 0.969238 | 0.0494719 | **+0.7347096** | 5.193637e−06 | −1.513742e−05 (291.5 %) | **−3.027485e−05** | **582.9 %** | −1.705500e−05 (328.4 %) |
| x = 3.80189e−04, Q² = 1.44699, y = 0.956277 | — | — | 4.824404e−06 | — | −1.294668e−05 | **268.4 %** | — |
| x = 3.46737e−04, Q² = 1.30202, y = 0.943487 | — | — | 4.485434e−06 | — | −6.667781e−06 | **148.7 %** | — |
| x = 3.16228e−04, Q² = 1.17158, y = 0.930869 | — | — | 4.174054e−06 | — | −4.106004e−06 | **98.4 %** | — |
| x = 2.88403e−04, Q² = 1.0542, y = 0.918419 | — | — | 3.828704e−06 | — | −2.759256e−06 | **72.1 %** | — |

Median over the 77 is **2.839e−05 %** of the band and the minimum
2.784e−10 %; **3** cells exceed 100 % of the band, **6** exceed 10 %, **8**
exceed 1 % and **18** exceed 0.1 %, so a handful of cells at the very top of
the y range carry essentially all of it. **Read the corner with
its health warning** — it is exactly where `TPeakPlusLL`'s own uncancelled soft
radiator `D(z) ∝ 1/(1−z)` makes the model untrustworthy (rc.hpp header block,
"two bad corners"), so a price of 583 % of the band there is not a licence to
quote it as a systematic; it is a statement that the omission is not small
where the model is not small.

**How much of §B1.3's collapse the bound recovers**, which is the only fair way
to score it:

| cell | `r_T/r_U` t-peak | t-peak+ll, sp = 0 | t-peak+ll, sp = 1 | **collapse recovered** |
|---|---|---|---|---|
| x = 0.01, Q² = 5 | −3.972572e−04 | −2.627428e−04 | −2.627428e−04 | **0.00 %** |
| x = 0.10, Q² = 5 | +4.161551e−05 | +1.324560e−07 | +1.324560e−07 | **0.00 %** |
| x = 0.30, Q² = 5 | +1.302390e−05 | +8.605653e−10 | +8.605653e−10 | **0.00 %** |
| x = 4.16869e−04, Q² = 1.6081 | −1.869598e−04 | −9.789112e−05 | −1.172405e−04 | **21.72 %** |
| x = 3.80189e−04, Q² = 1.44699 | −1.839451e−04 | −1.001863e−04 | −1.092468e−04 | 10.82 % |
| x = 3.46737e−04, Q² = 1.30202 | −1.810251e−04 | −9.998544e−05 | −1.050504e−04 | 6.25 % |
| x = 3.16228e−04, Q² = 1.17158 | −1.781957e−04 | −9.942787e−05 | −1.027674e−04 | 4.24 % |
| x = 2.88403e−04, Q² = 1.0542 | −1.754523e−04 | −9.786149e−05 | −1.002885e−04 | 3.13 % |

**21.72 % is the most this bound recovers anywhere on the production grid.**

**Whole run.** 200 000 events, seed 1234, config 1, `--rc-tail-model
t-peak+ll`, at the CLI's default fill **`P_z = 0.7`** (named 2026-09-15: the
`sp_tensor_scale = 0` row re-measures 1.040274188 / min 1.000002021 /
max 5.981595247 at `tensor_thirds_plan(0.7, 0.6)`, against 1.040368933 /
1.000002120 / 5.981595247 at `tensor_thirds_plan(0.0, 0.6)`):

| `sp_tensor_scale` | mean `rc_tail` | min | max | tail-clipped |
|---|---|---|---|---|
| **0 (default)** | **1.040274188** | 1.000002021 | 5.981595247 | 0 / 200 000 |
| 1 | 1.040274193 | 1.000002021 | 5.981588558 | 0 / 200 000 |
| 100 | 1.040274711 | 1.000002021 | 5.980926365 | 0 / 200 000 |

The mean barely moves — the term carries `Q_N/6` and the run averages over
m = ±1, 0, and it is alive on 77 of 3051 cells — but the **max** does, so it is
not an accounting fiction.

**Exactly linear**, so one run rescales (unlike `fq_scale`, which is quadratic
and must be run, T12): scale 0.5 / 2 / 4 give exactly 0.5× / 2× / 4× the
increment (`4.00000000000639` at the worst alive cell; rtol 1e−10 over the 6
production cells that clear the cancellation guard, worst 8.2435e−11), and the
increment is linear in `Q_N` too (−2 : +1 at 1e−9), which keeps T14's Eq. (43)
identity exact.

### B1.5 The ordering of the unpriced pieces, which is what the reader needs

At Q² = 5 and the three standard x, **as a fraction of the band half-width**:

| piece | x = 0.01 | x = 0.10 | x = 0.30 |
|---|---|---|---|
| the RC band itself (δ(x)·\|A_zz\|) | **100 %** by definition | 100 % | 100 % |
| `--rc-qe-tensor-scale 1` (polarised QE tail) | 0.0855 % | 0.00031 % | **0.246 %** |
| `--rc-sp-tensor-scale 1` (elastic s/p tensor) | **0.000 %** | **0.000 %** | **0.000 %** |
| the **quasi-elastic s/p** tensor column | **not priceable — no knob exists** | | |

At the low-x, `y → 0.97` corner the order **inverts and both blow past the
band**: `sp = 1` gives **582.9 %** of the band half-width against
`qe = 1`'s **328.4 %**.

**The gap this task leaves open, named.** `qe_tensor_scale` multiplies `σ^q_U`
alone; `sp_tensor_scale` multiplies `σ^el_{U,s+p}` alone. The **quasi-elastic
s/p column** `TailTriple::qe_sp` is therefore covered by **neither**, and it
keeps a tensor part of exactly zero. It is the column that *survives* where the
coherent one dies — by §B1.3, essentially **all** of the ×314 and ×15 134
denominator growth at x = 0.10 and 0.30 is quasi-elastic — so **the largest
exactly-zero tensor term in `rc_tail` is now the one no knob reaches.**

That is deliberate and not an oversight: covering it needs a **third** stand-in
that borrows across coherence *and* across Q'² in one product, i.e. both
category errors at once, and this task did not invent one. It is asserted as a
gap by `python/tests/test_rc.py::test_the_quasi_elastic_sp_column_is_bounded_by_NEITHER_scale`,
which shows that at the three standard points both scales at once leave the
tensor term exactly where `qe_tensor_scale` alone puts it, so the gap cannot be
quietly closed by accident or quietly forgotten.

### B1.6 The one departure from the brief, stated

The brief asked for the price *"at the standard points (x = 0.01, 0.10, 0.30;
Q² = 5) for scale 0, 0.5, 1, against the band half-width"*. That price is
**identically zero**, for the reason in §B1.4, so a table of three zeros is the
whole of the requested measurement and would have read as a knob that does
nothing. The requested table is published above **unaltered and unrounded**,
and a second table at the cells where the term is alive is published beside it,
labelled as such. Nothing was substituted for the requested number.

The brief also anticipated the failure mode — *"the s/p vertex sits at
Q'² ~ Q² where the coherent form factor is dead for ⁶Li, so the elastic tensor
fraction may be the WRONG reference entirely at high Q²"* — and asked that it
be confronted. Confronted: it is worse than a wrong reference. Where the
coherent form factor is dead the **term itself** is dead, so the knob does not
give a wrong price, it gives **no price**, and the piece that carries the tail
there is `qe_sp`, which this knob was never defined on.

### B1.7 What ships

| file | change |
|---|---|
| `include/lipolgen/rc.hpp` | `RcOptions::sp_tensor_scale` + its 135-line defensibility header ("A BOUND WITH NO DERIVATION", points (1)–(5)); the `RcTailModel::TPeakPlusLL` enumerator comment and the header block's tail-band bullet re-stated as "bounded, not computed" **with the numbers** |
| `src/core/rc.cpp` | the numerator term; the `>= 0` check; the `tail_model != TPeakPlusLL` refusal; the comment block at `tail_ratio_at` |
| `src/core/pipeline.cpp` | `validate()`'s `>= 0` check and refusal; `--rc-sp-tensor-scale` in the tail sub-knob reach list; the `knob_provenance` row with its three-way reach rule |
| `python/bindings.cpp` | `meta["rc_sp_tensor_scale"]`; the `RcOptions.sp_tensor_scale` docstring |
| `python/lipolgen/__init__.py` | `make_config(rc_sp_tensor_scale=…)`; `RC_TAIL_MODELS`' own "no tensor partner" note re-stated with the collapse factors (edited **line-neutrally** — two citations point into this file) |
| `python/lipolgen/cli.py` | `--rc-sp-tensor-scale`; the default; the pass-through; two run-banner branches carrying the measured "empty at the standard points / 583 % in the corner / qe_sp bounded by nothing" sentences |
| `tests/test_rc.cpp` | `TEST_CASE("B1: the s-/p-peak tensor stand-in…")` — 5 subcases: the default and the dead `TPeak` tables; the refusal pair and the negative; the numerator identity at 270 nodes (rtol 1e−5, measured worst 3.60359e−06); exact linearity in the scale and in `Q_N`; the price, with the empty standard points and the alive corner printed |
| `python/tests/test_rc.py` | 3 new tests + the `meta` key set, the `make_config` round trip and the CLI-switch list |
| `python/tests/test_knob_provenance.py` | 2 matrix cells (`read` at `t-peak+ll`, `refused` at `t-peak`) + the module docstring's refusal list |
| `tests/test_rc.cpp` (T8(d) subcase (d)) | its "UNKNOWN, not zero" comment re-stated as "BOUNDED … not computed", with the note that the bound is EMPTY at that subcase's own (x, Q²) |
| `docs/CONVENTIONS.md` | `rc_sp_tensor_scale` added to the REFUSE list of the "refuse or label" criterion, **line-neutrally** (§B1.8) |
| `validation/check_physics_channels_links.py`, `validation/physics_channels_ranges.json` | three allow-list line moves and the re-pointed citations/fingerprints of §B1.8 — **no allow-list entry added, none removed** |
| `docs/` | `SUMMARY.md` item 5, `OPEN_ITEMS_SOLUTIONS.md` §9 (two sites), `PHYSICS_CHANNELS.md` (a new row + "Ten" → "Eleven" CLI switches), `USAGE.md` §7b (a mandatory-band row + the prose, two sites) |

**No source file was added, so SPDX is unchanged at 102/102.** No physics
number is defined twice: the term reuses `ratio_t`, the same object the elastic
tensor term uses, and defines no new constant.

### B1.8 A bookkeeping defect found on the way

`docs/PHYSICS_CHANNELS.md`'s SND-configuration row cited
**`docs/USAGE.md:2152-2327` §9**. That range has not been §9 for some time — it
points into **§7b** (its recorded last line was *"the tensor part of those
peaks is **unknown, not zero**"*), while §9 *"Polarized ⁶Li configurations for
coherent-diffraction codes"* starts at line **2652** and runs to EOF at
**2849**. The fingerprint gate had
no way to see it, because it verifies the recorded first/last lines and those
still matched. It surfaced only because a B1 edit landed **inside** the range
and changed its last line. The citation is now **`docs/USAGE.md:2652-2849`**,
the section the row is actually about, and the range re-recorded.

**The rest of the line bookkeeping, stated in full.** B1's insertions grew
(`wc -l`) `rc.hpp` 1634 → 1788, `pipeline.cpp` 3639 → 3693, `rc.cpp`
1747 → 1781, `cli.py` 1532 → 1588 and `__init__.py` 836 → 840, which broke
**170** of the
gate's `file:line` citations. `check_physics_channels_links.py --fix` repaired
150 of them (the declaration-site ones). The remaining **18 distinct
citations** (19 occurrences in `PHYSICS_CHANNELS.md`, 14 keys in
`validation/physics_channels_ranges.json`) are **use sites**, which `--fix`
cannot resolve; each was re-pointed **by an old→new line map computed with
`difflib` against the pre-edit copy of its file**, so every one moved to the
line carrying the identical text, and the three allow-list entries inside the
checker itself (`pipeline.cpp:783 → 802`, `pipeline.cpp:1753 → 1772`,
`rc.cpp:1033 → 1053`) were moved by the same map.

`USAGE.md` §7b's own growth (its prose gained 6 lines and the mandatory-band
table 1) then shifted **two more range citations** — §8 *Command-line
generators* `2623-2630 → 2625-2632` and the freshly re-pointed §9
`2650-2847 → 2652-2849`. Both were moved **by content**: the citation and its
`physics_channels_ranges.json` key were renamed and the stored `sha256` left
untouched, so the gate re-verified the blocks rather than being told to accept
them — the same discipline §B3 used for the same two sections.

**And `docs/CONVENTIONS.md` had to be edited LINE-NEUTRALLY.** Its
"refuse or label" paragraph needed `rc_sp_tensor_scale` added to the REFUSE
list. `CONVENTIONS.md` is not itself a gated document, but eight range
citations in `PHYSICS_CHANNELS.md` and one in `SPIN32_FINITE_GAMMA.md` point
*into* it, and the first version of the edit added two lines and moved every
one of them (`66-74 → 68-76`, `194-198 → 196-200`, `372-374 → 374-376`,
`168-171 → 170-173`, `312-314 → 314-316`, `508-509 → 510-511`, and
`72-73` reported as "the last line of the range is blank"). Rather than
re-point nine citations for a three-clause sentence, the paragraph was
rewritten to fit the **same six lines** (`git diff --numstat` = `2 2`), which
leaves every one of them untouched. **The same line-neutrality was required
of two later edits and honoured**: `python/lipolgen/__init__.py`'s
`RC_TAIL_MODELS` note (two citations point into that file) and
`python/lipolgen/cli.py`'s `--rc-sp-tensor-scale` help text (three do), both
rewritten inside their existing line counts rather than re-pointing citations
for a comment. *(How this was caught, stated because the sequencing is the
lesson: the gate had been run and was green BEFORE the `CONVENTIONS.md` edit,
and the breakage surfaced in the next full `pytest` run, in
`python/tests/test_doc_link_gate.py`, which runs the same checker. The two
agree; what does not is a gate run and a later edit. Re-run the gate after
every docs edit, not once per phase.)*

`--record-ranges` **was** run here, twice and deliberately, unlike in §B3:
once for the re-pointed `docs/USAGE.md` §9 range above, and once to
fingerprint `src/core/rc.cpp:1010`, the new row's use-site citation of
`sp_tensor_scale` (the gate requires a fingerprint for a use site in a file
that names the symbol on more than one line). **No citation was deleted and no
cited block's claim was changed**; the reference count moved 1224 → **1228**,
which is exactly the four citations the new `PHYSICS_CHANNELS.md` row adds.

### B1.9 Suites and gates, before and after

Measured in this session, both ends. The "before" column is the working tree
as B1 found it — B3's edits present and uncommitted, so it is **not** §B0's
read-only baseline and is re-measured here rather than copied.

| check | before B1 | after B1 | delta |
|---|---|---|---|
| `build/lipolgen_tests` | **405 cases / 17 240 417 assertions / 0 failed / 1 skipped** | **406 / 17 241 655 / 0 failed / 1 skipped** | **+1 case, +1238 assertions** — the whole of `TEST_CASE("B1: …")` |
| `python -m pytest python/tests -q` | **1008 passed / 127 skipped** | **1016 passed / 136 skipped** | **+8 passed, +9 skipped** — 3 new `test_rc.py` tests plus the 14 new `test_knob_provenance.py` matrix cells (2 variants × 7 specs), of which **5 ran** and **9 skipped** because their base is refused on that spec's channel |
| `check_physics_channels_links.py` (strict) | **1224 refs / 96 ranges / 7 external / 0 broken / 6 allow-listed** | **1228 / 96 / 7 / 0 / 6** | **+4 refs** = the new `PHYSICS_CHANNELS.md` row's four citations |
| … `SPIN32_FINITE_GAMMA.md` | 19 / 115 / 0 broken | 19 / 115 / 0 broken | unchanged |
| … `PYTHIA_BRIDGE.md` | 0 / 8 external / 0 broken | 0 / 8 external / 0 broken | unchanged |
| `check_spdx_headers.py` | **102 files, all stamped** | **102 files, all stamped** | unchanged — **no source file was added** |
| `git diff --stat validation/reference/` | empty | **empty** | no rtol-1e-12 reference gate moved |

*(The docs-gate "before" number was re-measured on a **temporarily reverted**
tree: the first attempt ran while B1's source edits were already in place and
reported 170 broken citations from the line shift alone, which is a
measurement of the shift and not of the baseline. The four files were reverted
by inverse string replacement, the gate re-run to get 1224 / 0 broken, and the
edits restored from a byte-checked copy — md5 identical before and after.)*

### B1.10 Registry

**No new row.** B1 adds no author decision: `sp_tensor_scale` defaults to 0,
every default is bit-identical, and — like the six no-op opt-ins of
`../run_2026-09-03/AUTHOR_DECISIONS.md` §B15 (`c0_shape`, `tail_model`, `qe_tensor_scale`,
`a_transfer_frac`, `cluster_vmc_mc_sigma`, `pol_sf`) — it moves nothing at the
CLI default whichever way the author would go. It is the direct sibling of
`qe_tensor_scale`, which §B15 already lists as *not* a decision. If the author
wants the pair decided, §B15 becomes registry row 26; **nobody here made that
call**, and `STATUS.md` (this run's) records that no row was added.

**What B1 leaves open, for a registry that may want it later:** the
quasi-elastic s/p tensor column (§B1.5) is now the largest exactly-zero tensor
term in `rc_tail` and no knob reaches it.

---

## B2. `RcTailModel::PolradFull` — POLRAD Eq. (18) + Appendix B + Eq. (A.4), implemented and gated

**THE HEADLINE.** The exact tail exists, is opt-in, and **it is not inside the
band the two t-peak models were being run as.** Over the sampler's own 3051
accepted cells it lies between `t-peak` and `t-peak+ll` on **1725 (56.5 %)**
and outside on **1326 (43.5 %)**; the **ratio of the σ-weighted mean shifts**
is **0.428** (a ratio of means over the cells — *not* an event-weighted one,
and *not* a mean of the per-cell coverage fractions, which is 6.483 — see
§B2.5), and in the `Q² ≥ 20 GeV²`, `y ≤ 0.9` window it is **above both
edges**.
The two t-peak models are a **price range, not a confidence interval**, and
every site that said otherwise now says this.

**AND THE SECOND HEADLINE, which is a physics correction to a shipped
sentence.** `rc.hpp` said, in three places, *"POLRAD supplies no tensor s/p
peak and this file will not invent one."* That is true of **Eq. (38)** and
**false of the paper**: Eq. (18) + Eq. (A.4) carries the s- and p-peaks' own
tensor content, and `PolradFull` computes it. So `RcOptions::sp_tensor_scale`,
B1's bound-with-no-derivation, can now be **scored against a computed answer**
— §B2.5 — and it is refused on `PolradFull` for the **opposite** reason it is
refused on `TPeak`: not because the term did not run, but because it did.

**WHAT IS STILL NOT CLOSED, named in the same breath.** The **quasi-elastic**
tail is a sum over **spin-½ nucleons**, and Eq. (A.5) makes `Im_{5..8}`
identically zero: there is no tensor quasi-elastic structure function to put at
*any* of the three peaks. So the exact tail closes the **elastic** half of B1's
gap and leaves the quasi-elastic half exactly where B1 left it —
`qe_tensor_scale` is still the only stand-in, and it is still a borrowed
magnitude. `test_rc.cpp` T19(f) and `test_rc.py`'s new test assert the gap so
it cannot be quietly assumed closed.

### B2.1 What was built, and what Eq. (18) needs that Eq. (38) does not

`RcTailModel::PolradFull` threw from 2026-09-02 to 2026-09-06 (`rc.cpp`'s
constructor, `PipelineConfig::validate()`), and `tests/test_rc.cpp`'s **T9** was
the suite's one unconditional `doctest::skip(true)`. Both are gone.

New in `src/core/rc.cpp` (one anonymous-namespace block plus four public entry
points; **no new file**, so the SPDX count is unchanged at 102):

| piece | what it is |
|---|---|
| `polrad_tau_limits` | Eq. (14)'s exact `τ_{max,min} = (S_x ± √λ_Q)/(2M²)`, with the lower limit written `−Q²/(M² τ_max)` so no digit cancels |
| `full_kin` | the invariants **plus** the target polarisation four-vector: `η = 2(a_η k₁ + c_η p)`, `a_η = M/√λ_s`, `c_η = −S/(2M√λ_s)` — POLRAD `conkin`'s `+self,if=long` (adgh:779-781) — and `(ηq)`, `(ηK)` (adgh:788-789) |
| `full_level0` | Eq. (B.12)'s eight F-functions on Eq. (B.13)'s `B_{1,2}`, `C_{1,2}`, with **Eq. (B.14)** for `F_d` (τ = 0 is inside the range and the printed `τ^{-1}(C_2^{-1/2} − C_1^{-1/2})` is 0/0 there) |
| `full_lift1`, `full_lift2` | Eq. (B.8)'s two upper-index lifts, `F → F^η → F^{ηη}` |
| `full_ffu` | Eq. (B.4) + Eq. (B.5) (`T_{5j1} = T_{1j1}`, `T_{6j1} = T_{2j1}`) for the four rows this file needs. `T_3`, `T_4` are **not** transcribed: every term of them carries `P_L` and the whole A_zz programme runs at an unpolarised beam |
| `full_thetas` | Eq. (B.1)'s `θ_ij = Σ_k a_ik T_ijk` with Eq. (B.3)'s `a_ik` and Eq. (B.7)'s `q_ik` |
| `polrad_im_el_spin1` | **Eq. (A.4) itself**, split as `Im_i = u[i] + (Q_N/6) t[i]`. One definition, called by the quadrature, so **T9 gates the shipped contraction** |
| `tanh_sinh_quadrature` | double-exponential quadrature; the τ range is split at `τ_s = −Q²/S`, `0` and `τ_p = Q²/X` so each peak sits at a **panel edge**, where tanh-sinh clusters its nodes exponentially |
| `polrad_full_sigma_el` | the elastic pair `(σ_u^A, σ_q^A)` from **one** τ pass (they share `θ_ij`, which is 90 % of the cost) |
| `polrad_full_sigma_qe_u` | the same quadrature with Eq. (A.5)'s spin-½ `Im_{1,2}`, `Im_{5..8} ≡ 0`, and the same `pauli_suppression` |

**The five things Eq. (18) needs that Eq. (38) does not** are enumerated in the
header comment above `PolradFullPair`: the τ variable itself (which is what
brings the s-/p-peaks in), Appendix B's kernels, the polarisation four-vector
(Eq. (38) has no `η` at all — its tensor content is already contracted),
Eq. (A.4) at `Q_N ≠ 0`, and **the lepton mass**, which the t-peak forms drop
entirely and which is what regulates the peaks.

**Nothing is defined twice.** The nucleon form factors are `nucleon_ff`, the
spin-1 ones the same `Spin1ElasticFF`, the Pauli factor the same
`pauli_suppression`, and the `Q_N = 0` limit of `polrad_im_el_spin1` is gated
**against** `rosenbluth_spin1` rather than duplicating it.

### B2.2 FIVE transcription defects in polrad2t.tex's Appendix B — measured, not argued

`polrad_transcription_check.md` checked Eq. (18), Appendix **A** and Eq. (A.4)
line by line in 2026-09-02/03 and found the paper right except for `σ_u^C`'s
`Y₋`. **Appendix B had never been transcribed at all**, and it is wrong in five
places. Each was found the same way — the code did not reduce to Eq. (38) —
and each is settled by the FORTRAN (`adgh`, CPC `ADGH_v1_0`), which is
internally consistent.

> **CORRECTION, 2026-09-06 (verification pass) — two `adgh` LINE CITATIONS in
> this section pointed into the wrong deck.** No number and no transcription
> changes; only the line numbers were wrong, re-checked against
> `adgh_v1_0.gz` sha256 `83703668…1ea7b` (uncompressed sha256 `4a823438…11761`,
> 11 593 lines). `T_821`'s pairing: **`adgh:9155-9157` → `adgh:1161-1163`**
> (`tm3(6,2,n)` in `ffu`; 9155-9157 is `al2ll`'s own `write`/`end`).
> `4F_{2−}^{ηη}`: **`adgh:1082-1085` → `adgh:1085-1086`** (`eeir`). The same
> two, plus `adgh:9134-9137` → **`adgh:1134-1137`**, `adgh:1090-1091` →
> **`adgh:1093-1094`**, and the first-lift range `1075-1078, 1088, 1090` →
> **`1073-1075, 1089, 1092`**, were corrected in `src/core/rc.cpp` and
> `../run_2026-09-02/polrad_transcription_check.md` §10 on the same pass.

The measurement: with a Gaussian form factor **dead at the s-/p-peak vertex**
(so only the t-peak survives), the ratio `Eq. (18)_tensor / Eq. (38) σ_q^d` must
tend to 1 as `x_A → 0`. Deuteron, E = 27.6 GeV, y = 0.5, F_m sector, `b = 20
t_min`:

| reading | x_A = 0.003 | 0.006 | 0.012 |
|---|---|---|---|
| **`adgh` (shipped)** | **1.00313** | **1.00627** | **1.01257** |
| paper's `a_ik` (Eq. (B.3), no `1/M^{l_i−1}`) | −9345.05 | −4658.87 | −2315.75 |
| paper's `q_ik` (Eq. (B.7), `τ/M` inside `T_{ij2}`) | −167714 | −83651.6 | −41620 |
| paper's `T_821` pairing (Eq. (B.4), `S ηk₁ + X ηk₂`) | 13.9681 | 14.0043 | 14.0771 |
| paper's `4F_{2−}^{ηη}` (Eq. (B.8), no `τ`) | 7406.36 | 464.806 | 30.1189 |
| paper's `4F^{ηη}` (Eq. (B.8), no `s_η`) | −86047.5 | −42857.2 | −21261.8 |

1. **Eq. (B.3)'s `a_ik` are missing `1/M^{l_i−1}`.** The paper prints
   `{Q² − 3(ηq)², 6(ηq), −3}` for i = 5, 6 and `{ηq, −1}` for i = 4, 8; `adgh`
   has both divided by `M²` and `M`. Eq. (A.1) says why: the `Im_5` term is
   `g̃_{μν} k_n Im_5` with **`k_n = (3(qη)² − Q²)/M²`**, so `a_{51} = −k_n`, not
   `−M²k_n`. Confirmed independently by POLRAD's **Born** kernels
   (`bornin`, adgh:817-825): `tm(5) = −ek·tm(1)` with `ek = (3 apq² − Q²)/M²`
   exactly, and `tm(8) = (apq/M)·(…)`.
2. **Eq. (B.7)'s `q_ik` is applied at the wrong level.** The paper puts
   `q_{i2} = τ/M` *inside* `T_{ij2}`, so `a_{i2} = 6(ηq)` multiplies it too;
   `adgh` (adgh:1115-1116) adds a bare `τ/M²·T_{i,j−1,1}` **outside** `a_{i2}`.
   The paper's reading is wrong by five orders of magnitude.
3. **Eq. (B.4)'s `T_{821}` pairs `S` with `ηk₁` and `X` with `ηk₂` with a `+`.**
   The paper prints `(S ηk₁ + X ηk₂)F_{1+}`, i.e. `(ηK S_p + ηq S_x)/2`;
   `adgh`'s `tm3(6,2,n)` (adgh:1161-1163, in `ffu`) has `apn F_{1+} S_x + apq F_{1+} S_p`, i.e.
   `S ηk₁ − X ηk₂`. Note `T_{811}` two lines up pairs them the *other* way
   (`S ηk₂ + X ηk₁`) and **that** one the paper gets right, which is what makes
   this look like a sign slip rather than a convention.
4. **Eq. (B.8)'s `4F_{2−}^{ηη}` drops a `τ`.** The paper's own first-order
   relation is `2F_{2−}^η = (2F_d + F_{2+}) **τ** s_η + F_{2−} r_η`; the
   second-order line prints `(2F_d + F_{2+})(r_η s + s_η r)` with no `τ`.
   `adgh` (`tails`, adgh:1085-1086) keeps it.
5. **Eq. (B.8)'s `4F^{ηη}` drops an `s_η`.** The paper prints
   `F(r_η − τ s_η)² + 4F_i(r_η − τ s_η) + 4F_{ii}s_η²` — the middle term is
   **first** order in (r, s) where the other two are second. `adgh` has
   `4F_i(r_η − τ s_η)**s_η**`.

**None of these is a "the PDF is unreadable" artefact.** All five were read
from the arXiv **LaTeX source** `polrad2t.tex` (hep-ph/9706516), the same
source `polrad_transcription_check.md` §0 obtained; the fraction-ordering trap
that document warns about does not apply to any of them.

### B2.3 THE THREE GATES

**(i) The limit in which Eq. (18) MUST reduce to Eq. (38), stated.** Not "switch
the s-/p-peaks off by hand": Eq. (38) *is* POLRAD's ultrarelativistic t-peak
extraction of Eq. (18) (§2.1.3 B), so the two agree only where the other two
peaks are absent and only to `O(x_A)`. The limit is therefore **`x_A → 0` with
a form factor dead at the s-/p-peak vertex `Q'² ≈ z_s Q²`** — which isolates
the t-peak *without touching either quadrature*. Deuteron, E = 27.6 GeV,
y = 0.5, Gaussian `b = 20 t_min`, `S_A = A s`, `x_A = x/A`:

| sector | x_A | `Eq. (38) σ_u` | `Eq. (18) σ_u` | ratio | `Eq. (38) σ_q` | `Eq. (18) σ_q` | ratio |
|---|---|---|---|---|---|---|---|
| spin-0 (F_c only) | 0.003 | 6.483070e−04 | 6.514995e−04 | **1.00492** | −0.0 | −0.0 | *(identically 0)* |
| | 0.006 | 1.609831e−04 | 1.625761e−04 | **1.00990** | | | |
| | 0.012 | 3.970125e−05 | 4.049424e−05 | **1.01997** | | | |
| F_m only | 0.003 | 2.177393e−08 | 2.182404e−08 | **1.00230** | −1.980547e−08 | −1.986747e−08 | **1.00313** |
| | 0.006 | 2.169600e−08 | 2.179625e−08 | **1.00462** | −1.977532e−08 | −1.989933e−08 | **1.00627** |
| | 0.012 | 2.154051e−08 | 2.174097e−08 | **1.00931** | −1.971485e−08 | −1.996274e−08 | **1.01257** |
| F_q only | 0.003 | 2.615099e−11 | 2.632608e−11 | **1.00670** | 1.156366e−11 | 1.178013e−11 | **1.01872** |
| | 0.006 | 1.041630e−10 | 1.055656e−10 | **1.01347** | 4.474334e−11 | 4.646615e−11 | **1.03850** |
| | 0.012 | 4.131267e−10 | 4.243740e−10 | **1.02722** | 1.668991e−10 | 1.805316e−10 | **1.08168** |

The spin-0 column is `1 + 1.64 x_A` to three digits, which is the size an
ultrarelativistic `O(x_A)` error should be. **The three sectors are run
separately on purpose**: the full deuteron `σ_q^d` is a ~7× cancellation
between the F_m and F_q sectors, so *its* ratio is not a gate — it is the ratio
of a residual to a residual. `test_rc.cpp` T19(a) pins all nine numbers at
2e-4 **and** asserts the deviation shrinks monotonically as `x_A` falls
(1.15e-03 unpolarised / 1.57e-03 tensor at `x_A = 0.0015`).

**And the same gate fixes the NORMALISATION**, which was the one thing
`polrad_transcription_check.md` §8b could not settle from the paper: the ratio
tends to **1**, not to a constant, so Eq. (18) as written here is the
**whole-nucleus** `d²σ/(dx_A dy)` with no `1/A` of its own — exactly what
Eq. (38) is — and `RcModel` gives it the **same `1/A²`**. The paper's
`σ₁^el = (1/A) d²σ/dx_A dy` notation (polrad2t.tex:579) remains loose, and this
is now a measurement rather than an inference.

**(ii) T9, un-skipped.** `TEST_CASE("T9: the Q_N = 0 Rosenbluth limit of
Eq. (A.4), and the Q_N part")` replaces the `doctest::skip(true)` stub. It
asserts, through the **shipped** `polrad_im_el_spin1`: `Im_2|₀ = A(Q²)` and
`Im_1|₀ = B(Q²)/2` of `rosenbluth_spin1` at 1e-15 over seven `t`;
`Im_{3,4,5,6,7,8}|₀ == 0` **exactly** (they carry an explicit overall `Q_N`, so
anything weaker would let a tensor term leak into the unpolarised tail); the
six `Q_N` parts against a literal Eq. (A.4) written in the test file at 1e-14,
**with the corrected `Im_6`** (`4/(1+η_A)`, not `4η_A/(1+η_A)` — the defect
`polrad_transcription_check.md` §7.1 found in the design and fixed on
2026-09-04, now exercised for the first time); and that `Im_1`, `Im_2` are
**not** purely unpolarised, which is design §1.4.6's own warning.

**(iii) The transcription check's own tabulated numbers.** `T8(c')` already
gates the t-peak against `polrad_transcription_check.md` §8's three compiled
deuteron points. T19(c) extends it: at every one of them the exact tail must
be **positive** (T8(0)'s sign gate), must land on the **same side of zero** in
the tensor column (that column changes sign with x — the reason `σ^el_T` is a
separate table), and must **exceed** the t-peak, since it adds two peaks.
Measured full/t-peak, pinned at 2e-3:

| point | `σ_u` t-peak | polrad-full | ratio | `σ_q` t-peak | polrad-full | ratio |
|---|---|---|---|---|---|---|
| E = 27.6, x = 0.050, y = 0.60 | 4.153620e−05 | 4.422247e−05 | **1.0647** | 2.734422e−06 | 3.272269e−06 | **1.1967** |
| E = 27.6, x = 0.012, y = 0.50 | 1.966607e−03 | 2.312642e−03 | **1.1760** | 9.172112e−05 | 2.358360e−04 | **2.5712** |
| E = 11.0, x = 0.200, y = 0.50 | 1.222669e−07 | 1.891953e−07 | **1.5474** | −5.430094e−08 | −5.841761e−08 | **1.0758** |

**What is still NOT checked, said in the header, the banner, the CLI help and
the bindings.** The absolute normalisation is **not** validated against
Mo–Tsai or against any external exact tail: `[MT69]` is not in this tree and no
number from it is quoted anywhere. What is checked is **POLRAD-internal**
((i) and (ii)) and **against the leading-log fallback** (§B2.4). Nothing more.

### B2.4 Where PolradFull and TPeak(+LL) differ — measured

**The HERMES deuteron point (x = 0.012, y = 0.85, Q² = 0.5283)** — the point
four shipped documents quote as *"the t-peak is 23 % of the leading-log total,
low by 4.36×"*. The exact tail settles it:

| piece | `t-peak` | `t-peak+ll` | **`polrad-full`** | full/t | full/(t+ll) |
|---|---|---|---|---|---|
| elastic `σ^el_U` | 1.340571e−03 | 2.063992e−03 | **1.793760e−03** | 1.3381 | 0.8691 |
| quasi-elastic `σ^q_U` | 2.262505e−03 | 1.365507e−02 | **6.695847e−03** | 2.9595 | 0.4904 |
| **total** | 3.603076e−03 | 1.571906e−02 | **8.489606e−03** | **2.3562** | **0.5401** |
| elastic `σ^el_T` | 6.252324e−05 | *(no tensor s/p)* | **2.419056e−04** | 3.8691 | — |

**So the t-peak is low by 2.36×, not 4.36× — the leading-log edge overshoots by
1.85×**, and the overshoot is almost entirely the **quasi-elastic** peaks
(`ll_peaks_qe` is 2.04× the exact answer there while `ll_peaks_spin1` is
1.15×). This is the first number in this tree that prices the leading log
itself. T19(d) gates both ratios at 2e-3.

**⁶Li config 1, the standard points and the corners** (`x_A = x/6`,
`S_A = 6 s`, `s = 3980`, HO form factor, `n_eta = 128`):

| x | Q² | y | `σ^el_U` t-peak | polrad-full | ratio | `σ^el_T` t-peak | polrad-full | ratio | tensor fraction × |
|---|---|---|---|---|---|---|---|---|---|
| 0.01 | 5 | 0.1256 | 7.99822e−04 | 8.02845e−04 | 1.0038 | −2.44472e−06 | −2.48272e−06 | 1.0155 | **×1.0117** |
| 0.10 | 5 | 0.01256 | 4.79427e−07 | 4.99755e−07 | 1.0424 | +4.48610e−10 | +1.33621e−10 | 0.2979 | **×0.28574** |
| 0.30 | 5 | 0.004188 | 4.81229e−11 | 9.30535e−11 | 1.9337 | +3.83640e−12 | +6.11172e−12 | 1.5931 | **×0.82387** |
| 0.01 | 20 | 0.5025 | 9.93842e−04 | 9.94884e−04 | 1.0010 | −3.03776e−06 | −3.05076e−06 | 1.0043 | ×1.0032 |
| 0.10 | 100 | 0.2513 | 4.99606e−07 | 5.00649e−07 | 1.0021 | +4.68968e−10 | +4.55246e−10 | 0.9707 | ×0.96872 |
| 1e−3 | 1 | 0.2513 | 1.83245e−01 | 1.84894e−01 | 1.0090 | −2.75950e−04 | −2.97420e−04 | 1.0778 | ×1.0682 |
| **1e−4** | **0.2** | **0.5025** | 3.42247e+01 | 3.60491e+01 | 1.0533 | −3.32661e−02 | −8.28291e−02 | **2.4899** | **×2.3639** |
| **4.169e−4** | **1.608** | **0.9692** | 1.99171e+01 | 2.20808e+01 | 1.1086 | −2.48284e−02 | −7.64873e−02 | **3.0806** | **×2.7788** |
| 0.7 | 25 | 0.008973 | 1.26502e−18 | 1.29953e−18 | 1.0273 | −1.86730e−18 | −1.68555e−18 | 0.9027 | ×0.87869 |

Two things to read here. **The coherent s-/p-peaks are alive only at low Q'²**
— at Q² ≥ 5 the leading-log elastic column is *identically zero*
(`t+ll_u/tpk = 1.0000` at every row above except the two low-Q² ones), yet
`polrad-full` still finds +0.4 % to +93 % more elastic tail, because Eq. (18)'s
peaks carry `O(α)` non-log pieces and the exact kinematics that a single-`z`
collinear log does not. **And the tensor fraction goes UP, not down, in the
low-Q² corner** (×2.36, ×2.78): those are the cells where the coherent form
factor is still alive at `Q'²_s`, so the s-/p-peaks carry real tensor content
— exactly the term `sp_tensor_scale` bounds and cannot compute.

**The quasi-elastic column, ⁶Li at the shipped `k_F` = 0.169 GeV** (per
nucleon, `/A` applied):

| x | Q² | t-peak | t-peak+ll | polrad-full | full/t | full/(t+ll) |
|---|---|---|---|---|---|---|
| 0.01 | 5 | 6.27347e−06 | 2.08602e−05 | 1.40824e−05 | 2.2448 | 0.6751 |
| 0.10 | 5 | 3.65893e−08 | 1.56721e−05 | 7.39919e−06 | **202.22** | 0.4721 |
| 0.30 | 5 | 1.33701e−09 | 2.02526e−05 | 9.52409e−06 | **7123.4** | 0.4703 |
| 0.01 | 20 | 7.79526e−06 | 7.85505e−06 | 8.07188e−06 | 1.0355 | 1.0276 |
| 1e−3 | 1 | 6.83503e−04 | 1.95069e−02 | 9.78345e−03 | 14.314 | 0.5015 |
| 1e−4 | 0.2 | 8.26707e−02 | 2.30482e+00 | 1.07372e+00 | 12.988 | 0.4659 |

**The leading log is consistently ~2× the exact quasi-elastic answer** wherever
the peaks matter, and the t-peak alone is low by ×202 and ×7123 at x = 0.10 and
0.30 — and **every bit of that is tensor-blind on all three models.**

### B2.5 The per-cell and event-weighted census

Over the sampler's own **3051 accepted cells** (⁶Li config 1), on the
unpolarised ratio `r_U = w_tail − 1`, with `tail_max` raised to 1e30 so the
census measures the **model** and not the ceiling:

| comparison | cells > 1 % | % of σ | median \|Δ\| | worst | σ-weighted mean `w_tail` over the cells |
|---|---|---|---|---|---|
| `t-peak` → `t-peak+ll` | **268 / 3051 (8.78 %)** | 30.37 % | 3.013e−05 | **×1.594** at x = 4.16869e−04, Q² = 1.6081, y = 0.96924 | 1.02217080 → **1.04097066** (+1.8392 %) |
| `t-peak` → **`polrad-full`** | **112 / 3051 (3.67 %)** | 7.66 % | 4.814e−05 | **×4518.3** at x = 0.954993, Q² = 206.68, y = 0.054378 (6.648e−21 → 4.5173e+03) | 1.02217080 → **1.03021634** (+0.7871 %) |
| `t-peak+ll` → **`polrad-full`** | **103 / 3051 (3.38 %)** | 7.48 % | 3.544e−05 | ×4518.3, same cell | 1.04097066 → **1.03021634** (−1.0331 %) |

**Is it inside the band?** No:

* between the two edges on **1725 / 3051 (56.54 %)**, outside on **1326
  (43.46 %)**;
* the fraction of the `t-peak → t-peak+ll` gap it covers is **0.6555 median
  OVER THE 3027 CELLS WHOSE GAP IS NONZERO** (10 % 0.4697, 90 % 1367, min
  0.3577, max 2.674e+05). **Say the count:** on the other **24** cells
  `t-peak+ll` equals `t-peak` to the last bit, so the fraction is ±∞ there;
  over all **3051** the median is **0.6620**;
* **0.4280** is the **RATIO of the σ-weighted mean shifts**,
  `(⟨full⟩−⟨t-peak⟩)/(⟨t-peak+ll⟩−⟨t-peak⟩)` on the σ-weighted means of the
  last column above — **0.427958** with the ceiling raised, **0.427368** at the
  shipped `tail_max = 10`. It is a **ratio of means, not a mean of the per-cell
  fractions**: the σ-weighted mean of those 3027 fractions is **6.483**, 15×
  larger, because a handful of tiny-gap cells carry fractions up to 2.7e+05;
* per **event** (200 000, seed 1234) it is outside on **10 046 (5.02 %)** at
  the CLI's default fill **P_z = 0.7**, and on **10 043 (5.02 %)** at
  **P_z = 0** — an event count is a statement about a fill plan.

*(Every number in the bullets above re-measured 2026-09-15 on the
working-tree build, at both P_z; the cell census itself is P_z-free — the
sampler's cell grid and `cell_xsec_pb` are identical on the two plans,
measured. Two published digits moved with the re-measurement: the 90 %
quantile, printed here as **1368**, is **1366.95**; and the median was printed
without the count it is taken over, which is the defect this bullet now
closes.)*

*(The last column of that table weights each accepted CELL by its own
`cell_xsec_pb`; the whole-run table below weights actual EVENTS. They are two
estimators of the same thing and they land 5.1e−04 apart — 1.03021634 against
1.029702912, i.e. 4.99e−04 relative (re-measured 2026-09-06; the "4e-4" that
stood here is neither the absolute nor the relative gap) — which is the size of
the sampler's own cell-to-event reweighting, not a disagreement. **Say the P_z
on the event side** (added 2026-09-15): 1.029702912 is the default fill
`P_z = 0.7`; against the `P_z = 0` run's 1.029775347 the same cell estimator
lands **4.41e−04** apart, **4.28e−04** relative. The cell estimator itself does
not move with P_z.)*

> **CORRECTION RECORD, 2026-09-06 (verification pass) — the shipped headers
> were quoting a THIRD number, mislabelled.** `include/lipolgen/rc.hpp:265` and
> the `RcTailModel::PolradFull` enum comment, `python/lipolgen/cli.py`'s
> `--rc-tail-model` help and `tests/test_rc.cpp`'s T19(g) comment carried
> *"Event-weighted mean rc_tail 1.02217080 → 1.04097066 → 1.03020526"* and
> *"×11 the t-peak at x = 0.955, Q² = 186"*. Neither is what this section
> publishes, and both were re-measured here before being corrected at every
> site:
>
> * **1.02217080 / 1.04097066 / 1.03020526 are CELL-σ-weighted, not
>   event-weighted, and the third is CLIPPED.** Re-measured over the same 3051
>   accepted cells: cell-σ-weighted **1.022170796 / 1.040970661 /
>   1.030205260** at the shipped `tail_max = 10`, and **1.030216345**
>   unclipped — the figure the table above publishes. The **event**-weighted
>   triple (200 000, seed 1234) is **1.021778527 / 1.040274188 / 1.029702912**
>   **at the CLI's default fill P_z = 0.7**, and **1.021836305 / 1.040368933 /
>   1.029775347** at **P_z = 0** (added 2026-09-15: this bullet named no P_z,
>   and the two differ in the fourth decimal place); it is identical to the
>   whole-run table below, which is the P_z = 0.7 run. The cell-σ-weighted
>   triple is the same on both plans.
> * **×11 is `1 + tail_max`, the CEILING, not the model.** Re-measured: at the
>   shipped `tail_max = 10` the worst `polrad-full`/`t-peak` ratio is exactly
>   ×11 and **six** cells tie at it — all at x = 0.954993, Q² = 185.977,
>   206.683, 229.696, 389.396, 432.752, 733.631, together **4.96e−09** of the
>   cross section and **0** of 200 000 events. Unclipped the model's own worst
>   is **×4518.26** at x = 0.954993, **Q² = 206.683**, y = 0.0544; the
>   Q² = 186 cell the headers named is **×3449.85** there.

**The worst cell is real physics and a numerical nothing at once.** At
x = 0.955 the DIS Born is negligible while the elastic/quasi-elastic tail is
not, so `w_tail` is genuinely ~4500 there; the t-peak gives **6.6e−21**,
because `t_min ∝ x²` has killed both the coherent and the nucleon form factor
at its own vertex while the s-/p-peaks sit at far lower `Q'²`. It is **not** a
lower bound with a 10 % error there — it is zero where the answer is
everything. That cell carries **5e−9 of the cross section** and **0 of 200 000
events**, and at the shipped `tail_max = 10` it is clipped:

| model | `clipped_cell_fraction` (table nodes) | by y-band (`<0.5` / `0.5-0.9` / `>0.9`) | accepted cells at the ceiling | events clipped (200 k) |
|---|---|---|---|---|
| `t-peak` | 0 | 0 / 0 / 0 | 0 / 3051 | 0 |
| `t-peak+ll` | 0.00360036 | 0 / 0 / **0.046205** | 0 / 3051 | 0 |
| **`polrad-full`** | **0.00385753** | **0.0045576** / **0.0012376** / 0 | **6 / 3051** (4.96e−07 % of σ) | **0** |

Note the y-bands: `t-peak+ll` clips at `y → 1` (its s-peak beating `Y₊`) and
`polrad-full` clips at `y → 0` — a **different corner**, and the same one where
`TPeakPlusLL`'s uncancelled soft radiator breaks down.

**Whole run, 200 000 events, seed 1234, ⁶Li config 1 inclusive, AT THE CLI's
DEFAULT FILL `P_z = 0.7` (`tensor_thirds_plan(0.7, 0.6)`):**

| model | mean `rc_tail` | min | max | tail-clipped | wall clock |
|---|---|---|---|---|---|
| `t-peak` (default) | **1.021778527** | 1.000000035 | 3.410153365 | 0 / 200 000 | 2.1 s |
| `t-peak+ll` | **1.040274188** | 1.000002021 | 5.981595247 | 0 / 200 000 | 2.2 s |
| **`polrad-full`** | **1.029702912** | 1.000001677 | 4.295059618 | 0 / 200 000 | **11.8 s** |

**The same run at `P_z = 0`** (`tensor_thirds_plan(0.0, 0.6)` — the plan the
C++ and pytest suites use, and the plan the 2026-09-03 run published), same
build, same seed, re-measured 2026-09-15:

| model | mean `rc_tail` | min | max |
|---|---|---|---|
| `t-peak` (default) | **1.021836305** | 1.000000035 | 3.410153365 |
| `t-peak+ll` | **1.040368933** | 1.000002120 | 5.981595247 |
| **`polrad-full`** | **1.029775347** | 1.000001650 | 4.295059618 |

The fill plan does not touch `RcModel` — it changes which events the sampler
draws, so **every event-weighted mean in this section is a statement about a
P_z as well as a seed**. The cell-σ-weighted census above is P_z-free
(measured: identical cell grid, `cell_xsec_pb` and `tail_ratio_at` on both
plans).

`x`, `q2` and `weight` are **byte-identical** across all three: `RcModel::fill`
takes no `Rng&`, so no tail model can move the random stream (T6's rule).

**In the `Q² ≥ 20 GeV²`, `y ≤ 0.9` window**, both fills, re-measured
2026-09-15 (5194 of 200 000 events at `P_z = 0.7`; **5182** at `P_z = 0`):

| model | `<w_tail − 1>`, `P_z = 0.7` | vs t-peak | `<w_tail − 1>`, `P_z = 0` | vs t-peak |
|---|---|---|---|---|
| `t-peak` | 8.719649e−03 | — | 8.489138e−03 | — |
| `t-peak+ll` | 8.773752e−03 | +0.6205 % | 8.540716e−03 | +0.6076 % |
| **`polrad-full`** | **8.815031e−03** | **+1.0939 %** | **8.581236e−03** | **+1.0849 %** |

**The exact tail is ABOVE both edges in that window**, so even the
event-weighted statement the band was quoted for does not bracket it.

> **CORRECTION, 2026-09-15 — THE 2026-09-06 NOTE THAT STOOD HERE WAS A FALSE
> ATTRIBUTION, AND IT IS WITHDRAWN.** That note read: *"A stale number
> corrected in passing. `rc.hpp`'s header block, `USAGE.md` §7b and
> `PHYSICS_CHANNELS.md` all quoted this window as "5182 of 200 000 events,
> 8.48914e−03 → 8.54072e−03". Re-measured, it is 5194 events and
> 8.719649e−03 → 8.773752e−03. The move is not B2's: the pre-B2 build's own
> 200 k npz (`w_tpeak.npz`, generated during B1) already gives
> 5194 / 8.719649e−03, so it happened earlier in this run — Phase A moved the
> ⁶Li kernel normalisation, which moves the sampler's own cell weights and
> therefore which events land in the window."*
>
> **Nothing moved.** The 2026-09-06 re-measurement ran the CLI/pipeline at its
> **default fill `P_z = 0.7`**; the "5182 / 8.48914e−03 → 8.54072e−03" it
> declared unreproducible is the **`P_z = 0`** row, which is the fill those
> numbers were published on. Re-measured on the working-tree build at both
> (the table above): **5182** and **8.489138e−03 → 8.540716e−03** at
> `tensor_thirds_plan(0.0, 0.6)` — the 2026-09-03 digits, exactly — against
> **5194** and **8.719649e−03 → 8.773752e−03** at
> `tensor_thirds_plan(0.7, 0.6)`. `w_tpeak.npz` gives 5194 for the same reason:
> re-read 2026-09-15, it carries mean **1.021778527**, window **5194**,
> **8.719649e−03** — the P_z = 0.7 run to every digit, so it too was generated
> at the default fill, and it was never evidence that anything had moved.
>
> **Phase A is exonerated.** It moved neither this window nor the run mean; the
> 2026-09-03 run's own record had already tabulated both rows
> (`run_2026-09-03/phase_B_numbers.md:836-839`, *"Why `--pz 0` is in that
> command"*), which is where this should have been checked before a cause was
> published. The correction, with both fills named, is carried at every site
> that quoted the note: `rc.hpp`'s header block, `USAGE.md` §7b and its
> `--rc-tail-model` row, `PHYSICS_CHANNELS.md`, `OPEN_ITEMS_SOLUTIONS.md` §9
> and §B1's own paragraph above.

**The tensor fraction `r_T/r_U`, and B1's bound scored against a computed
answer.** At the three standard points, Q² = 5 (pipeline-level `tail_ratio_at`,
so the quasi-elastic term and the `Q_N/6` factor are in it):

| x | `t-peak` | `t-peak+ll` | **`polrad-full`** | t+ll × | **full ×** |
|---|---|---|---|---|---|
| 0.01 | −3.972572e−04 | −2.627428e−04 | **−3.154008e−04** | ×0.66139 | **×0.79395** |
| 0.10 | +4.161551e−05 | +1.324560e−07 | **+8.336376e−08** | ×0.0031829 | **×0.0020032** |
| 0.30 | +1.302390e−05 | +8.605653e−10 | **+2.942224e−09** | ×6.6076e−05 | **×0.00022591** |

**B1's `sp_tensor_scale` was not even one-sided.** At x = 0.01 and 0.30 it
pointed the right way and stopped far short (it recovered **0.0 %** of the
collapse at both, being bit-identical to 0 there, while the computed answer
recovers **39.1 %** and **×3.4** of `t-peak+ll`'s residue); at x = 0.10 the
computed tensor fraction is **below** `t-peak+ll`'s, i.e. the bound pointed the
**wrong way**. Over all 3051 cells the `polrad-full`/`t-peak` tensor-fraction
multiplier runs **−0.098877 (min) / 0.002266 (10 %) / 0.90719 (median) /
0.9975 (90 %) / 348.53 (max)**, and it is **higher** than `t-peak+ll`'s on
**1886 / 3051 cells (61.8 %)**. Read `sp_tensor_scale` as a price tag on an
omission — which is what its own field comment says — and not as an interval.

**The price on `ΔA_zz`.** `Δ(ΔA_zz) = 2[r_T(model) − r_T(t-peak)]/(1 + r_U)`,
band half-width `δ(x)·|A_zz|` on this tree's own ⁶Li b₁, Q² = 5:

| x | band half-width | `t-peak+ll` − `t-peak` | % of band | **`polrad-full`** − `t-peak` | **% of band** |
|---|---|---|---|---|---|
| 0.01 | 1.358472e−04 | 1.200916e−07 | 0.0884 % | **5.908189e−08** | **0.0435 %** |
| 0.10 | 9.680016e−05 | 6.316227e−07 | 0.6525 % | **2.974583e−07** | **0.3073 %** |
| 0.30 | 7.307492e−07 | 5.053477e−08 | 6.915 % | **2.377717e−08** | **3.254 %** |

At the three standard points the exact tail's price is **about half** the
leading-log edge's, and the band remains ×768–×1086 wider than either (B6).

### B2.6 Numerics: two things that are refused rather than degraded

**(a) `long double`, and it is not a style choice.** Eq. (B.3)'s `a_ik` for
i = 5, 6 are the expansion of `−M²k_n(q − k) = Q² − 3(η(q−k))²` in powers of the
photon contraction, and in the collinear region the three terms cancel:
measured `Σ|term| / |Σ term|` inside one `θ_5j` is **2.9e4** at ⁶Li config 1,
x = 1e−4, y = 0.5 near `τ_min`, and **6.6e4** at x = 0.1 — ~4.8 decimal digits
gone before the i-sum starts. In plain `double` the **tensor** column then
moves by **14 %** between two tanh-sinh step sizes the unpolarised column
agrees on to ten digits (−3.30e−04 at h = 1/16 against −3.78e−04 at 1/32,
x = 1e−3, y = 0.5): the quadrature is converged and the **arithmetic** is not.
In `long double` (64-bit mantissa) the same scan is stable to **6e−05**.
`RcModel`'s constructor **refuses** `PolradFull` where
`numeric_limits<long double>::digits <= numeric_limits<double>::digits`,
rather than ship a tensor tail with two significant figures; T19 asserts the
platform property.

**(b) `n_eta ≥ 64`, because it means something different here.** Under
`PolradFull` `n_eta` is the **tanh-sinh node count per panel** of the τ
integral (four panels, split at `τ_s`, 0, `τ_p`). Measured at ⁶Li config 1,
x = 1e−3, y = 0.5:

| `n_eta` | `σ^el_U` | `σ^el_T` |
|---|---|---|
| 32 | 0.2197427678 | −3.447295e−04 |
| 64 | 0.2208599758 | −3.447606e−04 |
| **128 (default)** | **0.2208659734** | **−3.446895e−04** |
| 256 | 0.2208659734 | −3.447235e−04 |
| 512 | 0.2208659734 | −3.447291e−04 |

`σ^el_U` is **0.51 % low at 32** and 2.7e−05 low at 64, against a value stable
to 1e−12 from 128 up; `σ^el_T` is stable to **~6e−04**, which is what the
cancellation of (a) leaves. `RcModel` and `PipelineConfig::validate()` both
refuse `n_eta < 64` on this model. T19(e) gates 1e−8 on `σ^el_U` and 1e−3 on
`σ^el_T` between 128 and 512.

**Cost.** The tail-table build (101 × 77 nodes, four columns) goes from
**0.15 s** to **9.84 s**; a 200 k-event run from **2.1 s** to **11.8 s**. Opt-in.

### B2.7 `m_lepton` — a reservation discharged

`RcOptions::m_lepton` was `Refused` unconditionally, with the reason *"RESERVED
for tail_model = PolradFull (POLRAD's F_IR and l_m). Nothing in
src/core/rc.cpp reads it"*. `PolradFull` reads it now — it is the `m²` of
Eq. (B.13)'s `C_{1,2}(τ)`, of `F_IR = m²F_{2+} − Q_m²F_d`, of `Q_m²` and of
`λ_s`, i.e. **what sets the width of the s- and p-peaks** — so the row is
three-way: `read` under `polrad-full`, `refused` on the two t-peak models,
`not-read` at `--rc off`. Measured, e → μ at ⁶Li x = 1e−4, y = 0.9:
`σ^el_U` **161.318 → 129.527 (×0.802932)** — a heavier lepton narrows the peaks
and the tail falls.

> **A second, smaller defect fixed on the way.** At `--rc off` the row said
> `refused`, and that was **false**: the whole rc validation block sits inside
> `if (rc != PipelineRc::Off)`, so `validate()` **accepts** any `m_lepton`
> there (verified). It now says `not-read` with the same sentence every other
> rc knob uses at `--rc off`. This is the one `meta` change on an `--rc off`
> run; the numeric columns are untouched (§B2.8).

### B2.8 The default did not move — measured, not argued

The `t-peak` path was not restructured (the `PolradFull` branch is an early
return placed **before** the existing lines), but bit-identity was still
checked empirically, against **the pre-B2 build's own artefacts**:
`w_tpeak.npz` and `w_0.0.npz`, the 200 000-event seed-1234 npz B1 generated
from the pre-B2 library.

| run (200 k, seed 1234, ⁶Li config 1 inclusive) | numeric columns | `meta` |
|---|---|---|
| `--rc tensor-band` (`t-peak`, the default) | **51 / 51 byte-identical** | no key added or removed, no non-provenance value changed, `knob_provenance` **69 → 69** rows, **2 rows changed and both only in their `reason` text** (`rc_m_lepton`, `rc_tail_model`); **no status and no value moved** |
| `--rc tensor-band --rc-tail-model t-peak+ll` | **51 / 51 byte-identical** | identical to the row above |
| `--rc off` (2000 events, seed 4242, against a 2026-09-04 artefact) | **48 / 48 byte-identical** | `rc_m_lepton` moves `refused` → `not-read` (§B2.7's defect fix); everything else unchanged |

*(The `--rc off` comparison is against an artefact predating B1 **and** B3, so
its 48/48 is a stronger statement than B2 alone needs — nothing in this run's
three phases moved an `--rc off` column.)*

`git diff --stat validation/reference/` is **empty**: no rtol-1e-12 reference
gate moved.

### B2.9 Suites and gates, before and after

The "before" column is the working tree as B2 found it — B3's and B1's edits
present and uncommitted.

| check | before B2 | after B2 | delta |
|---|---|---|---|
| `build/lipolgen_tests` | **406 cases / 17 241 655 assertions / 0 failed / 1 skipped** | **408 / 17 241 949 / 0 failed / 0 skipped** | **+2 cases** (T19 and the rewritten T9), **+294 assertions**, and **the suite's last `doctest::skip` is gone** |
| `python -m pytest python/tests -q` | **1016 passed / 136 skipped** | **1023 passed / 151 skipped** | **+7 passed, +15 skipped** — 1 new `test_rc.py` test plus 21 new `test_knob_provenance.py` matrix cells (3 variants × 7 specs), of which 6 ran and 15 skipped because their base is refused on that spec's channel |
| `check_physics_channels_links.py` (strict) | 1228 / 96 / 7 / **0 broken** / 6 allow-listed | **1236** / 96 / 7 / **0 broken** / 6 allow-listed | **+8 citations** on the `PHYSICS_CHANNELS.md` row that owns the tail (the four new entry points, header and source); ranges, externals and allow-list unchanged |
| … `SPIN32_FINITE_GAMMA.md` | 19 / 115 / 0 broken | 19 / 115 / 0 broken | unchanged |
| … `PYTHIA_BRIDGE.md` | 0 / 8 external / 0 broken | 0 / 8 external / 0 broken | unchanged |
| `check_spdx_headers.py` | 102 files, all stamped | **102 files, all stamped** | unchanged — **no source file was added** |
| `git diff --stat validation/reference/` | empty | **empty** | no rtol-1e-12 reference gate moved |

**The line shifts were absorbed, and one of them needed hands.** `rc.cpp` grew
by ~470 lines and `pipeline.cpp`, `cli.py`, `__init__.py` and `rc.hpp` by less;
`--fix` re-pointed **192** of the 213 shifted citations automatically. The
remaining **21** are **use-site pins** — citations to a line that does not
declare its symbol — which `--fix` cannot relocate. Each was re-pointed **by
content**: the exact line text is recorded in
`validation/physics_channels_ranges.json`, so every one was found by searching
for its own recorded line and the fingerprint re-verified rather than
re-blessed. Three of them are `ALLOW`-table entries in
`check_physics_channels_links.py` itself, whose keys carry the line number, and
those three keys were moved with them (`pipeline.cpp:802 → 841`,
`pipeline.cpp:1772 → 1811`, `rc.cpp:1053 → 1521`). **No citation was deleted
and no cited block was edited.**

> *(2026-09-06, verification pass: the `rc.cpp` allow-table key moved once more,
> `1521 → 1532`, and eight use-site pins were re-pointed by content and
> re-recorded — `rc.cpp:1424/1470/1495/1514 → 1435/1481/1506/1525`,
> `__init__.py:412/624/684 → 416/628/688`, `cli.py:1411 → 1427` — with every
> recorded line text and sha256 **unchanged** (`--record-ranges` reported
> `8 new, 0 changed, 8 dropped`). See §B2.12.)*

### B2.10 What B2 does NOT do — the list, so nobody has to infer it

1. **It is not validated against Mo–Tsai**, or any external exact tail. Said in
   the enumerator comment, the run banner, the CLI help, the bindings and here.
2. **It does not close the quasi-elastic tensor gap.** Eq. (A.5) has
   `Im_{5..8} ≡ 0` for a spin-½ target, so the quasi-elastic tail is
   tensor-blind at all three peaks on all three models — and by §B2.4 it is
   **202×** and **7123×** the t-peak at x = 0.10 and 0.30, i.e. it is where the
   tensor fraction's collapse actually lives. `qe_tensor_scale` remains the
   only stand-in and remains borrowed. **Gated** by T19(f) and by
   `test_rc.py::test_polrad_full_is_the_exact_tail_and_does_not_bracket_the_t_peak_band`.
3. **It is not the default**, and no published number moves (§B2.8).
4. **It carries no vector (P_N) sector.** `T_{3jk}`, `T_{4jk}` and
   `Im^el_{3,4}` are not transcribed: every term of them reaches the tail
   through `m M P_L` and the whole A_zz programme runs at an unpolarised beam,
   which `RcModel` already refuses to violate.
5. **It is longitudinal-only in `η`.** The axis dependence is carried, exactly
   as under `TPeak`, by `Q_N → P_zz^eff = 3 Q_NN P₂(cos θ_S)` (Eq. (43)); the
   tensor sector is **quadratic** in `η`, so its sign does not enter, but a
   genuinely transverse `η` (POLRAD's `+self,if=tran`) is **not** transcribed.
6. **It does not touch the second-order `alpha2ll` corrections**, the
   inelastic radiative tail (Eq. (13)), or the virtual/soft `δ` that POLRAD's
   own total carries — `rc_tail` is and remains the **hard-photon tail alone**.

### B2.11 Registry

**No new row.** `PolradFull` defaults to nothing — `TPeak` is still the shipped
tail and every published number is byte-identical — so, like the six no-op
opt-ins of `../run_2026-09-03/AUTHOR_DECISIONS.md` §B15 (`c0_shape`,
`tail_model`, `qe_tensor_scale`, `a_transfer_frac`, `cluster_vmc_mc_sigma`,
`pol_sf`), it moves nothing at the CLI default whichever way the author would
go. It is a **new value of the `tail_model` knob §B15 already lists**, not a
new decision. `../run_2026-09-06/STATUS.md` records that no row was added.

**What B2 leaves for a registry that may want it later, and it is now sharper
than B1 could put it:** the choice of tail model is no longer "a lower bound
and a stated model" — there is a computed answer, it is **outside** the pair on
43.5 % of cells, and the default is still the lowest of the three. Whether the
shipped default should stay `t-peak` when `polrad-full` exists is an author
decision; **nobody here made it**, and the price of making it is measured in
§B2.5 (event-weighted `rc_tail` 1.021778527 → 1.029702912 at the default fill
P_z = 0.7 and 1.021836305 → 1.029775347 at P_z = 0 — **+0.776 %** and
**+0.777 %**, re-measured 2026-09-15; the "+0.79 %" this line printed is the
**cell-σ-weighted** pair's +0.7871 % from the table above, not the
event-weighted one's).

### B2.12 The 2026-09-06 verification pass — what two adversarial reads found, and what was done

Two lenses were run over B1/B2/B3 after they landed. **B3's constant, B1's
bound header, and every term of the Eq. (18)/Appendix B/Eq. (A.4)
transcription survived** — sixteen single-term perturbations of `src/core/rc.cpp`
were rebuilt and run, and fourteen of them BITE. Four residues did not, and all
four are closed here.

**1. `RcTailModel::PolradFull` had no pinned VALUE anywhere.** T19(a)–(h) and
T9 pin `polrad_full_sigma_el` / `polrad_full_sigma_qe_u` directly; T14 and
`test_rc_pipeline` check only refusals; the only assertion that touched the
model's own branch of `RcModel::tail_sigma_at` was T19(g)'s inequality
`tt.u > mt.tail_sigma_at(0.01, 5.0).u`. **Measured**: replacing that branch's
per-nucleon `1/A²` by `1/A` — the *exact* historical factor-A bug, at the exact
line it lived on (§4 of `../run_2026-09-02/polrad_transcription_check.md`) —
left **T9 + T19 at 316/316** and both `polrad-full` pytests green.

Closed by **T19(i)** (`tests/test_rc.cpp`) and a mirror block in
`python/tests/test_rc.py::test_polrad_full_is_the_exact_tail_and_does_not_bracket_the_t_peak_band`,
which pin `RcModel::tail_sigma_at` to VALUES at three points on ⁶Li config 1
(`RcOptions()` defaults but for the model; σ **per nucleon**, GeV⁻²):

| x | Q² [GeV²] | y | `u` = σ^el_U | `t` = σ^el_T | `qe` = σ^q_U |
|---|---|---|---|---|---|
| 0.01 | 5 | 0.12563 | **2.2301252856028388e−05** | **−6.8964416598016638e−08** | **1.4082400494855565e−05** |
| 0.10 | 20 | 0.050251 | **1.3474277075174561e−08** | **1.0398154385161343e−11** | **4.2214155843320782e−08** |
| 0.30 | 5 | 0.0041876 | **2.5848200248542323e−12** | **1.6976986824052815e−13** | **9.5240888233400158e−06** |

pinned at rtol **1e−9**, plus the same three points against the raw
quadratures at rtol 1e−12 (`u·A² == polrad_full_sigma_el().u`,
`qe·A == polrad_full_sigma_qe_u()`). The first row is the reviewer's point:
`2.2301252856028388e−05` is `8.0284510281702199e−04 / 36`, **not** `/6`. The
second is where `x_A = x/A` and `x` sit 7.8 decades apart on the form factor;
the third is where `qe` is **3.7e+06 ×** `u`, so the single 1/A of the
quasi-elastic column is pinned where it is the whole tail.

**The pin bites, measured both ways** (each perturbation rebuilt, run, and
`src/core/rc.cpp` restored to sha256 `24a2238e…19dea`):

| perturbation in the `PolradFull` branch | T9 + T19 before T19(i) | T9 + T19 with T19(i) | polrad-full pytests |
|---|---|---|---|
| `per_nucleon` `1/A²` → `1/A` | 316/316, **0 failed** | **12 failed** (u, t and both ratio pins × 3 points) | **1 failed** (`0.00013380751713617034` against `2.2301252856028388e−05`) |
| `x_a` → `x` at the elastic call | 1 failed (T19(g) only) | **13 failed** | — |

**2. B2's own physics correction was not propagated.** *"POLRAD supplies no
tensor s/p peak"*, which §B2 itself declares true of Eq. (38) and **false of
the paper**, still stood unqualified at eleven sites. All are now qualified
("Eq. (38) supplies none; Eq. (18) + Eq. (A.4) does, and `PolradFull` computes
it"): `rc.hpp` (`sp_tensor_scale`'s *WHY IT EXISTS*, and `TailTriple`, which
was a **third** occurrence of the *"unknown, not zero"* phrase B1 reported
closed), `bindings.cpp`, `cli.py`'s `--rc-sp-tensor-scale` help (which
contradicted the `--rc-tail-model` help twenty lines above it),
`__init__.py`'s `RC_TAIL_MODELS` note, `src/core/rc.cpp`,
`PHYSICS_CHANNELS.md`'s row (line 807 of the same file already said the
opposite), `tests/test_rc.cpp`'s B1 case, `python/tests/test_rc.py`'s B1
docstring, and §B1 of this document. The three **dated** records
(`../run_2026-09-02/polrad_transcription_check.md` §8,
`../run_2026-09-02/phase_C_numbers.md`, `../run_2026-09-03/phase_B_numbers.md`
§B2.1) got correction boxes rather than edits. `design_C_tensor_rc.md` §3.1/§3.2
still specified `PolradFull = 1 … NOT IMPLEMENTED` and the matching
`validate()` refusal; both now carry dated notes.

**3. The published census was the ceiling, mislabelled.** See the correction
record in §B2.5: the shipped triple was cell-σ-weighted and clipped, printed as
"event-weighted", and ×11 is `1 + tail_max`. Corrected at
`rc.hpp:265`, the `RcTailModel::PolradFull` enum comment, `cli.py`'s
`--rc-tail-model` help and T19(g)'s comment.

**4. Transcription nits, each re-verified against the FORTRAN** (`adgh_v1_0.gz`
sha256 `83703668…1ea7b`; uncompressed `4a823438…11761`, 11 593 lines):

| site | was | is |
|---|---|---|
| `rc.cpp` `FullFSet` header, `ffu` rebuilding `F_IR` | `adgh:9134-9137` (`al2ll`'s `refss=0d0`) | **`adgh:1134-1137`** (`hi2`/`shi2`/`ehi2`/`ohi2` in `ffu`, `adgh:1125`) |
| `rc.cpp` + §B2.2 + `polrad_transcription_check.md` §10.4, `T_821`'s pairing | `adgh:9155-9157` (`al2ll`'s `write`/`end`) | **`adgh:1161-1163`** (`tm3(6,2,n)`) |
| `rc.cpp` `full_lift1` | `adgh` `tails`:`1075-1078, 1088, 1090` (1076-1078 is `ois`, the η–s CROSS set) | **`1073-1075, 1089, 1092`** (`eis`, `eir`, `ei12`; `ei1pi2`; `eb`) |
| `rc.cpp` `full_lift2` | `adgh` `tails`:`1082-1093` | **`1083-1094`** (`eeis` … `eeb`) |
| `rc.cpp` `full_level0` | `adgh` `tails`:`1053-1069` | **`1053-1065`** (1066-1069 is `sps`/`spe`/`ccpe`/`ccps`, i.e. s_η and r_η) |
| §B2.2 + §10.5, `4F_{2−}^{ηη}` keeps its `τ` | `adgh:1082-1085` | **`adgh:1085-1086`** (`eeir`) |
| §10.6, `4F^{ηη}` keeps its `s_η` | `adgh:1090-1091` | **`adgh:1093-1094`** (`eeb`) |
| `rc.hpp` spin-0 Eq. (18)/Eq. (38) ratios | 1.00493 / 1.00990 / 1.01998 | **1.00492 / 1.00990 / 1.01997** — re-measured 1.0049243949, 1.0098956443, 1.0199740696 at n_τ = 128 **and** 512, which is what T19(a) and §B2.3 already carried |
| `rc.hpp` Appendix-B defect count | "wrong in **three** places" | **FIVE** — §B2.2 and §10.2–10.6, one subsection per defect |

**End state, re-measured on this pass** (⁶Li config 1 unless stated):

| check | after B2 | after the verification pass |
|---|---|---|
| `build/lipolgen_tests` | 408 cases / 17 241 949 assertions / 0 failed / 0 skipped | **408 / 17 241 975 / 0 failed / 0 skipped** (+26 = T19(i)) |
| `python -m pytest python/tests -q` | 1023 passed / 151 skipped | **1023 passed / 151 skipped** (the pin is a block inside an existing test) |
| `check_physics_channels_links.py` | 1236 / 96 / 7 / 0 broken / 6 allowed | **1236 / 96 / 7 / 0 broken / 6 allowed** — 162 citations re-pointed (152 by `--fix`, 10 by content), 1 allow-table key moved, 8 pins re-recorded |
| `check_spdx_headers.py` | 102 | **102** |
| `git diff --stat validation/reference/` | empty | **empty** |
| 200 k, seed 1234 mean `rc_tail` (fill `P_z = 0.7`, the CLI default — named 2026-09-15) | 1.021778527 / 1.040274188 / 1.029702912 | **identical**, and the `--rc off` 48-column run is unchanged |

**No default moved and no code changed.** Every edit outside `tests/` and
`python/tests/` is a comment, a docstring or a `help=` string; `src/core/rc.cpp`
differs from its pre-pass state **only in comment lines** (checked line by
line), and the shipped tail model is still `TPeak`.

---

## B3. Q(⁷Li) sourced, and the A = 7 quadrupole gate made a real, cited test

`SUMMARY.md` item 4: Q(⁷Li) = −4.00(3) fm² was **quoted from memory and absent
from this tree**. It is now a sourced constant with one home, and the **A = 7
quadrupole gate** of `../run_2026-09-03/phase_D_li7_rank2.md` §4.2 is a
committed test that
computes Q from the in-tree α–t overlap and **reports a ratio**.

### B3.1 What was fetched, and what was read in it

Network was used for this item only. Three documents, all published and free:

| URL | HTTP | what was read |
|---|---|---|
| `https://nucldata.tunl.duke.edu/nucldata/ourpubs/07_2002.pdf` | 200, 573 975 B | TUNL's A = 7 half of Tilley, Cheves, Godwin, Hale, Hofmann, Kelley, Sheu, Weller, *"Energy Levels of Light Nuclei, A = 5, 6, 7"*, **Nucl. Phys. A708 (2002) 3**. The ⁷Li GENERAL block prints, verbatim: `µ = +3.256427(2) nm: see (1989RA17).` and **`Q = −40.6 ± 0.8 mb (1988DI1B).`**, followed by *"See (1988DI1B) for a review of earlier determinations, particularly those of (1984SU09, 1984VE03, 1984VE08, 1985WE08)"*. Its own reference list gives **1988DI1B = G.H.F. Diercksen, A.J. Sadlej, D. Sundholm, P. Pyykkö, Chem. Phys. Lett. 143 (1988) 163**, **1985WE08 = A. Weller *et al.*, Phys. Rev. Lett. 55 (1985) 480** (the "Weller" of the task brief), **1984SU09 = D. Sundholm, P. Pyykkö, L. Laaksonen, A. Sadlej, Chem. Phys. Lett. 112 (1984) 1** |
| `https://nucldata.tunl.duke.edu/nucldata/ourpubs/06_2002.pdf` | 200, 1 015 692 B | The **A = 6 half of the same paper**, read only to confirm the convention match: it prints `µ = +0.8220473(6) nm, +0.8220567(3) nm` and **`Q = −0.818(17) mb (1998CE04)`** — i.e. exactly what `rc.hpp`'s existing comment on `LI6_QUADRUPOLE_FM2` says the TUNL A = 6 evaluation prints, digit for digit |
| `https://www.nndc.bnl.gov/nndc/stone_moments/nuclear-moments.pdf` | 200, 555 008 B | N. J. Stone, *"Table of Nuclear Magnetic Dipole and Electric Quadrupole Moments"* (NNDC mirror, 152 pp., dated 04/11/2001). Page 1 carries the whole ⁷Li block (eight Q entries, listed in `rc.hpp`) and the ⁶Li entry `−0.00083(8) st [7Li] MB,R CPL 112 1 (84)`; the *Policies* section supplied the two rules that make the block readable — Q is **in barns**, and *"where two values of Q are given based on CER experiments, the first represents the value assuming constructive interference … the second assumes destructive"* |

`https://www-nds.iaea.org/publications/indc/indc-nds-0833.pdf` (Stone's 2021
IAEA recommended-value table) returned **HTTP 402** and was not read; nothing
is quoted from it, or from Stone 2016 / Pyykkö, which were not fetched.

### B3.2 The value, and why this one

```
LI7_QUADRUPOLE_FM2 = -4.06        (include/lipolgen/rc.hpp, beside LI6_QUADRUPOLE_FM2)
```

−40.6 ± 0.8 mb × 0.1 fm²/mb = **−4.06(8) fm²**, spectroscopic convention
(the m = I substate; negative = oblate).

**The choice is provenance, not physics.** The brief says to match the
convention of the tree's own Q(⁶Li) = −0.0818 fm², which came from the TUNL
A = 6 evaluation. TUNL serves A = 5, 6 and 7 as three PDFs of **one paper**,
so taking the A = 7 half gives both lithium quadrupoles from a single
document, a single sign convention and a single unit rule. The registry's own
D11 row had already framed the choice this way (`−4.00(3) / −4.06 fm², with a
source`, `LI7_QUADRUPOLE_FM2 beside LI6_QUADRUPOLE_FM2 in rc.hpp`) and had
already priced it at `ratio 0.871 → 0.858`; that price is confirmed below.

Stone lists **eight** ⁷Li determinations and recommends none (all in barns,
transcribed into `rc.hpp`): −0.0406 st (MB,R, CPL 112 1 (84)), −0.0370(8)
(CIAN, PRL 55 480 (85) — Weller), −0.041(6) (OD,OL), −0.059(8) (OL),
−0.040(11) (CER), **−0.0400(6)** and **−0.0400(3)** (CER, NP A530 475 (91) —
Voelk, the constructive- and destructive-interference readings of one
experiment), −0.0406(8) (R, AuJP 42 597 (89)). The note's −4.00(3) fm² is the
seventh of those. The spread is −3.70 … −5.9 fm²; the two values in play here
differ by **1.5 %**.

### B3.3 What was built

| where | what |
|---|---|
| `include/lipolgen/rc.hpp:669` | `LI7_QUADRUPOLE_FM2 = -4.06`, with the provenance above, the eight-entry Stone table, the swap price, and the sentence that it is a **gate reference and not a model input** |
| `include/lipolgen/b1_nuclear.hpp` | `struct AlphaTQuadrupole` (`r2_fm2`, `r_rms_fm`, `z_eff`, `q_fm2`), `alpha_t_quadrupole(const VmcRadial&, …)`, `li7_alpha_t_quadrupole(vmc_mc_sigma, r_max_fm, n_r, n_k)` |
| `src/core/b1_nuclear.cpp` | the Fourier–Bessel L = 1 transform, Z_eff from `nuclear_mass`, and a **j₁ branch added to the file-local `sph_bessel`** — the l = 0 and l = 2 branches are untouched, so every number this file already published is bit for bit what it was |
| `python/bindings.cpp` | `LI7_QUADRUPOLE_FM2`, `AlphaTQuadrupole`, `alpha_t_quadrupole`, `li7_alpha_t_quadrupole` |
| `tests/test_b1_nuclear.cpp` | **T13**, the gate: 26 assertions |
| `python/tests/test_li7_rank2.py` | **G8**, four tests |

plus the documentation sites (§B3.6). **Every source edit is an ADDITION**:
`git diff --numstat` reads `+74/−0` `b1_nuclear.hpp`, `+53/−0` `rc.hpp`,
`+24/−0` `bindings.cpp`, `+106/−0` `test_b1_nuclear.cpp` and `+62/−1`
`b1_nuclear.cpp` — the single deleted line being `sph_bessel`'s own error
string, `"... implemented for l = 0, 2"` → `"... for l = 0, 1, 2"`.

**No file under `src/` other than `b1_nuclear.cpp` changed at all**, so
`rc.cpp` — the whole radiative-correction path — is byte-identical and
`--rc off` cannot have moved; `git diff --stat validation/reference/` is
**empty**. No option, no CLI flag, no `knob_provenance` row: nothing here is a
knob, so there is no reach rule to state and no matrix cell to add. The new
constant is read by exactly one function, and that function is called only by
tests.

### B3.4 The measurement

All of it from `data/vmc/momenta/li7_at3.momentum` (the 3/2⁻ **ground** state)
through `li7_alpha_channel(BETA_DEFAULT, VmcAV18, σ)`, at the shipped defaults
r_max = 30 fm, n_r = 4000, n_k = 8001.

| quantity | measured 2026-09-06 |
|---|---|
| ⟨r²⟩ | **12.534828961030 fm²** |
| r_rms | **3.540456038568 fm** |
| Z_eff = Z_α(M_t/M₇)² + Z_t(M_α/M₇)² | **0.695075101362** (the A-number form 34/49 = 0.693878 is **0.172291 %** away and is *not* used) |
| Q = −(2/5) Z_eff ⟨r²⟩ | **−3.485059004257 fm²** |

**It reproduces the offline note to the last bit.** Running
`phase_D_li7_rank2.md` §9's own numpy `q_at` in the same interpreter against
the shipped C++ gives **relative difference 0.000e+00** on all three of
⟨r²⟩, Q and Z_eff — the same double, not "agreement to N digits". (The C++
uses the same natural cubic spline through `ClusterPartialWave`, the same
trapezoid and the same grids; the j₁ series threshold differs and does not
reach the printed digits.)

Grid/truncation, and the ANL Monte Carlo band (fully correlated,
`vmc_mc_sigma` = ∓1), both **pinned in T13 and G8**:

| r_max [fm] | Q [fm²] | | σ | Q [fm²] | vs σ = 0 |
|---|---|---|---|---|---|
| 20 | −3.472004271671 | | −1 | −3.496426569222 | **+0.326180 %** |
| **30 (default)** | **−3.485059004257** | | 0 | −3.485059004257 | — |
| 40 | −3.485797071382 | | +1 | −3.473744615578 | **−0.324654 %** |
| 60 | −3.486171274781 | | | | |

30 → 60 fm moves Q by **0.031915 %**; n_r/n_k 4000/8001 → 8000/16001 moves it
by **4e−9** relative (−3.485059004 → −3.485059017). So neither the r cut nor
the quadrature is anywhere near the discrepancy below.

### B3.5 The gate: a ratio, reported

| against | ratio | low by |
|---|---|---|
| `LI7_QUADRUPOLE_FM2` = **−4.06** (TUNL, adopted) | **0.858388917305** | **14.1611 %** |
| −4.00 (Voelk CER, Stone's seventh entry) | **0.871264751064** | **12.8735 %** |

The note's own **13 %** is the second row; the adopted constant gives **14 %**.
The compilation choice moves the ratio by **1.50 %** and does **not** change
the verdict, which is *"the α–t wave function reproduces the measured
quadrupole to 13–14 %, converged to 0.03 %, with a ±0.33 % MC band"*. Both
ratios are pinned, at both sites.

**The contrast is the result, and it is measured too.** The same construction
one cluster level down —
`alpha_d_quadrupole_fm2(li6_alpha_d_partial_waves())` = **−0.333315 fm²**
against `LI6_QUADRUPOLE_FM2` = −0.0818 — misses by a factor **4.0748**, i.e.
**307.476 %** against ⁷Li's **14.1611 %**: a factor **21.7**. "The ⁷Li
wave-function input is validated by a measured moment an order of magnitude
better than ⁶Li's is" is therefore a measurement, not a manner of speaking,
and T13/G8 assert it with a threshold of 10.

### B3.6 What this does NOT claim — said at every site

Written into `rc.hpp`'s constant comment, `b1_nuclear.hpp`'s declaration
comment, T13's header comment, T13's own `MESSAGE`, the `AlphaTQuadrupole` and
`li7_alpha_t_quadrupole` binding docstrings, G8's docstrings, and every
document site that carried the old flag — `docs/PHYSICS_CHANNELS.md`'s
spin-3/2 tensor-slot row, `docs/USAGE.md` §2c, `docs/OPEN_ITEMS_SOLUTIONS.md`
row 15 / §15.4 / §15.5 D11, `docs/benchmarking/03_data_nuclear.md` §5,
`run_2026-09-03/`'s `SUMMARY.md` item 4, `phase_D_li7_rank2.md`
(§1 headline, §4.2, §7 D11, §8 step 2, §10), `phase_D_numbers.md` and
`AUTHOR_DECISIONS.md` §B20:

* it validates the **α–t WAVE FUNCTION's quadrupole** — its ⟨r²⟩ and its
  P-wave character — and nothing else;
* it says nothing about the light-cone convolution or the DIS input, which is
  where the real ⁷Li uncertainty lives;
* **b₁(⁷Li) remains unimplemented** (open item 15; it waits on registry row 3 /
  D2, the unpolarised-backend decision) and ⁷Li's whole rank-2 sector is still
  exactly zero by construction. `test_the_passing_wave_function_gate_ships_no_b1_at_all`
  runs the gate and then re-asserts that zero in the same test, so a passing
  wave-function gate can never be read as licensing a b₁;
* there is **no A = 3 analogue** of the A = 2 b₁ gate and there cannot be one
  (³H and ³He are J = ½), so this quadrupole is the only offline validation the
  ⁷Li wave function will get;
* S_αt = 1.0084 > 1 proves the α–t overlap is not a probability; the missing
  13–14 % is the expected size of cluster distortion and non-α–t components.

### B3.7 Suites and gates, before and after

| check | B0 baseline (`53f6951`) | after B3 |
|---|---|---|
| `build/lipolgen_tests` | 404 cases / 17 240 391 assertions / 0 failed / 1 skipped | **405 cases / 17 240 417 assertions / 0 failed / 1 skipped** (+1 case = T13, +26 assertions) |
| `python -m pytest python/tests -q` | 1004 passed / 127 skipped | **1008 passed / 127 skipped** (+4 = G8) |
| `check_physics_channels_links.py` | 1220 / 96 / 7 / 0 broken / 6 allowed; SPIN32 19/115/0; PYTHIA 0/0/8 | **1224** / 96 / 7 / **0 broken** / 6 allowed; SPIN32 19/115/0; PYTHIA 0/0/8 — +4 citations, everything else identical |
| `check_spdx_headers.py` | 102 files, all carry it | **102 files, all carry it** (no new source file; the new code went into existing files) |

**The +4 is exactly what this task cited**: `LI7_QUADRUPOLE_FM2`,
`AlphaTQuadrupole`, `alpha_t_quadrupole` and `li7_alpha_t_quadrupole`, added to
`PHYSICS_CHANNELS.md`'s spin-3/2 tensor-slot row (the row that already said
⁷Li's rank-2 sector is identically zero, which is where the gate belongs and
where it is now stated). Ranges, externals, broken and allow-listed are
unchanged.

The line shifts were absorbed rather than papered over. `rc.hpp` grew 53 lines
and `python/bindings.cpp` 24, moving **60** named citations;
`check_physics_channels_links.py --fix` re-pointed every one — **58** in
`docs/PHYSICS_CHANNELS.md` (55 into `rc.hpp`, 3 into `bindings.cpp`), spread
over 13 table rows, and **2** range citations in
`docs/theory/SPIN32_FINITE_GAMMA.md` with their 2 fingerprint records. Editing
`docs/USAGE.md` §2c then moved **6** of that file's own range citations: 4
`--fix` re-pointed, and **2 were re-pointed by CONTENT by hand** — §7 *"Checks
you get for free"* and §8 *"Command-line generators"*, which `--fix` reported
as *"the first/last line of the range is blank"* rather than as moved. For
those two the citation and its `validation/physics_channels_ranges.json` key
were moved to the block whose recorded `first` and `last` lines they still
match exactly (`2008-2014` → `2014-2020`, `2581-2588` → `2587-2594`, both a
+6 shift), and **the stored `sha256` was left untouched**, so the gate
re-verified the content rather than being told to accept it. **No citation was
deleted, no cited block was edited, no allow-list entry was added, and
`--record-ranges` was never run** — it would have re-blessed drift instead of
reporting it.

### B3.8 Registry

**No new row.** This pays the existing **D11** of `AUTHOR_DECISIONS.md` §B20
(`STATUS.md` decision row 20 / `OPEN_ITEMS_SOLUTIONS.md` §15.5), whose
"what changes when decided" column already named `LI7_QUADRUPOLE_FM2` beside
`LI6_QUADRUPOLE_FM2` in `rc.hpp` and §4.2's gate, and whose predicted price
(`ratio 0.871 → 0.858`) is confirmed to six digits. `../run_2026-09-06/STATUS.md`
records that no row was added.
