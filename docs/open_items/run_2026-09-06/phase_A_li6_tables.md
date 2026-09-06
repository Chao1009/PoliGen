# Run 2026-09-06 — Phase A4: the ⁶Li b₁ tables, regenerated as a BAND

Baseline this section started from, re-measured before anything moved (tree
carrying A1–A3): **404 doctest cases / 17 240 384 assertions / 1 skipped /
0 failed**; **992 pytest passed / 127 skipped**; docs gate **1220 references
(strict) / 95 ranges / 7 external / 0 broken / 6 allow-listed**, +19/115
(`SPIN32_FINITE_GAMMA.md`), +8 external (`PYTHIA_BRIDGE.md`); SPDX **101/101**.

**What this section does, and what it does not decide.**
`run_2026-09-02/phase_D_numbers.md` carries the ⁶Li four-term convolution
tables; every row of them was made with the **default `ToyF2`** unpolarised
input, and `OPEN_ITEMS_SOLUTIONS.md` §10 has said since 2026-09-03 that they
*"must be REGENERATED under `--b1-unpol mstw` before they are quoted as
physics"*. **They are regenerated here — not on one backend, but as the BAND
over the two selectors that reach them.** No default moved, no option was
added, and no number in `validation/reference/*.json` can move: this section
constructs `Li6ConvolutionB1` objects and runs `Pipeline`s, and writes nothing.
`--b1-unpol` and `--unpol-sf` both already existed and are both already
registered in `Pipeline::knob_provenance`.

> ## The three conditions that ride with EVERY number below
>
> **(i) The mandatory ±100 % band.** Every ⁶Li b₁ value here is the **centre**
> of a `{0, 1, 2} × b₁` band and may only be quoted as that band. This does
> **not** come from the A = 2 gate and is **not** lifted by anything in this
> document: it comes from **Q(⁶Li) = −0.0818(17) fm²** against
> **Q_d = +0.2859(3) fm²** — the α–d D wave enters the closest measured
> observable with the **opposite sign to the deuteron's** and nearly cancels it.
> **Never quote one row.**
>
> **(ii) The A = 2 gate says nothing about the α–d step.** `DeuteronConvolutionB1`
> validates the *kernel* on the deuteron. There is **no measurement of b₁ at
> any A > 2**, so the whole α–d construction — terms (2d), (2α) and (3) — is
> ungated. Condition (i) exists because of condition (ii).
>
> **(iii) Quote the configuration.** Unpolarised backend **first**, then
> `r_sigma_lt` and κ = 1 (`finite_q_delta = false`), which are the ⁶Li
> backend's defaults and **not** the gate's (`r1998`, CDKS Eq. (21)).

---

## A4.1 The recipe — like for like with the tables being replaced

`phase_D_numbers.md`'s tables were made from **`Li6ConvolutionB1` directly**
(that is what `validation/b1_li6_table.py` does: `four_terms()` builds a bare
`lg.Li6ConvolutionB1()` and evaluates it at Q² = 2.5 on
x = 0.05/0.10/0.20/0.30/0.50). Everything below is the **same object at the
same points with the same options**, with one field set:

```python
import lipolgen as lg
lg.lhapdf_quiet()
BACKENDS = [("toy", None),                    # --b1-unpol toy   (the default)
            ("ct18nlo", lg.LhapdfSF("CT18NLO")),   # --b1-unpol ct18nlo
            ("mstw", lg.MstwSF())]                 # --b1-unpol mstw

def model(unpol, **kw):
    o = lg.Li6ConvolutionOptions()
    if unpol is not None:
        o.unpol = unpol                       # <- the ONLY field this varies
    for k, v in kw.items():
        setattr(o, k, v)
    return lg.Li6ConvolutionB1(o)

Q2 = 2.5
for name, u in BACKENDS:
    b = model(u)
    for x in (0.05, 0.10, 0.20, 0.30, 0.50):
        print(name, x,
              x * b.b1_embedded_s(x, Q2),          # term (1)
              x * b.b1_cg_dwave(x, Q2),            # term (3)
              x * b.b1_alpha_d_dwave_d(x, Q2),     # term (2d)
              x * b.b1_alpha_d_dwave_alpha(x, Q2), # term (2a)
              x * b.b1(x, Q2, 0.0))                # total
```

`o.unpol = None` (the default) resolves to `ToyF2` inside the backend, and an
explicitly constructed `lg.ToyF2()` is **bit-identical** to it — measured, so
the "toy" column below is literally `phase_D_numbers.md`'s own column and not
a re-derivation of it. **Every `toy` value in this document reproduces
`phase_D_numbers.md` to every printed digit**, which is what makes the other
two columns a like-for-like band and not a different measurement.

The second axis is the pipeline's, and is run through the CLI's own
constructor:

```python
cfg = lg.make_config(isotope="6Li", channel="inclusive", config=1,
                     events=2000, seed=99, b1_model="li6-convolution",
                     b1_unpol=B1U, unpol_sf=USF)     # B1U, USF in {toy, ct18nlo, mstw}
g = cfg.grid; g.nx, g.nq2 = 12, 6; cfg.grid = g      # SPEED only
sc = cfg.scenario; sc.x_max = 0.95; cfg.scenario = sc
p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
p.dis_sampler.kernel.tables(x, 2.5)                  # .b1, .b2, .f1, .f2
```

**`sc.x_max = 0.95` and the coarse grid interact, and the claim is stated with
the grid it was measured on** (an earlier draft of this section said the
refusal fires "at the shipped `x_max = 1.0`" without naming a grid, which is
true on the shipped grid and false on the one this recipe runs).

* On the **shipped grid** (`nx = 100`, `nq2 = 72`) at the shipped
  `x_max = 1.0`, **all seven legal cells — the `toy`/`toy` default included —
  are refused by `InclusiveSampler`**: *"negative phi-averaged density for
  m=1 at x = 0.955, Q2 = 167.3"*, with 1 + w_avg = **−0.8916** (`toy`/`toy`),
  −0.8373 / −12.71 / −11.56 (`ct18nlo` under `toy`/`ct18nlo`/`mstw`) and
  −0.842 / −12.75 / −11.59 (`mstw` under the same three). That is the
  pre-existing refused-base condition the pytest matrix already skips 127
  cells for, and `--x-max 0.95` is what `USAGE.md` §2a prescribes for this
  backend. (The other two of the nine cells never reach the sampler: they are
  refused earlier, by `validate()`, for the provenance reason below.)
* On the **coarse 12 × 6 grid this recipe uses for speed**, x = 0.955 is not a
  cell centre, so nothing is refused at either `x_max` and the flag is
  **inert**: σ is **equal as doubles** at `x_max` 0.95 and 1.0 in all seven
  legal cells (measured; e.g. `toy`/`toy` 5.98929947653798736e+05 both ways).
  It is carried here so the configuration matches the shipped guidance, not
  because it changes a number in this document.

Either way it is a property of the sampler window and has nothing to do with
either selector.

This document was produced by scripts run from `env.sh`; they build no state
and are reproduced verbatim above. Both optional tiers must be present
(`lg.HAVE_LHAPDF` and `lg.HAVE_PYTHIA8` both `True`, with the CT18NLO set and
`mstw2008lo.00.dat` on disk); a build without one **cannot** reproduce the
matching column and `validate()` refuses the selector rather than downgrade it.

---

## A4.2 Which axis moves what — MEASURED, not assumed

The band has two axes and they are **not** two views of one thing.

| axis | what it is | what it moves in these tables |
|---|---|---|
| `--b1-unpol` (`Li6ConvolutionOptions::unpol`) | the **deuteron F₁ inside CDKS Eq. (22)**, i.e. the unpolarised function the α–d *orbital* terms are folded against | **terms (2d) and (2α), and nothing else** |
| `--unpol-sf` (`InclusiveKernel::Options::f2_source`) | the **kernel's own F₁/F₂**, i.e. the **denominator** D_φ of A_zz and the unpolarised rate | **nothing at all in these tables**; it moves A_zz and σ |

Both statements are measurements.

**The `--unpol-sf` axis is exactly flat on x·b₁, and here is why it must be.**
`Li6ConvolutionB1::b1(x, q2, f1)` **ignores its `f1` argument** — the parameter
is unnamed in the definition (`src/core/b1_nuclear.cpp:710`) — so the kernel's
own F₁ has no path into b₁ except through `Li6ConvolutionOptions::unpol`, and
`default_inclusive_kernel` fills that slot with the caller's `b1_unpol` object
whenever one is named (`o.unpol = b1_unpol ? std::move(b1_unpol) : f2;`,
`src/core/pipeline.cpp:1168`). The one case where the kernel's object *does*
leak in is `b1_unpol` **null** — and that is precisely the combination
`PipelineConfig::validate()` refuses. Measured over the whole matrix: for
`b1_unpol` = `ct18nlo` and for `mstw`, the kernel's `tables(x, 2.5).b1` at all
five x is **bit-identical as doubles** across `unpol_sf` ∈ {toy, ct18nlo,
mstw}, and each equals the direct `Li6ConvolutionB1` value **bit for bit**.

**The 3 × 3 matrix: 7 cells run, 2 are refused by name.**

| `--b1-unpol` \ `--unpol-sf` | `toy` | `ct18nlo` | `mstw` |
|---|---|---|---|
| **`toy`** | runs (the shipped default) | **REFUSED** | **REFUSED** |
| **`ct18nlo`** | runs | runs | runs |
| **`mstw`** | runs | runs | runs |

The refusal is `PipelineConfig::validate()`'s, by name, and it is a
**provenance** refusal rather than a physics one:

> `PipelineConfig: unpol_sf = mstw with b1_unpol = toy on b1_model =
> li6-convolution: `toy` there means "the kernel's own UnpolSF, shared as one
> object", and that object is no longer ToyF2 -- so meta["b1_unpol"] would
> record "toy" for a b1 folded against mstw.`

So the two refused cells are **not** missing measurements: they are the two
cells in which `meta["b1_unpol"]` would lie. The band over x·b₁ is therefore
**three-wide, not seven-wide**, and that is a fact about the code, established
above, not an economy.

**What `--unpol-sf` does move**, at x = 0.10, Q² = 2.5 (same runs):

| `--unpol-sf` | kernel F₁ | kernel F₂ | σ (pb, 2000 ev, seed 99, `x_max` 0.95) |
|---|---|---|---|
| `toy` | 1.489486725 | 0.3489654613 | 5.989299477e+05 |
| `ct18nlo` | 1.662154064 | 0.3894189522 | 4.911979215e+05 (×0.820126) |
| `mstw` | 1.481358591 | 0.3470611556 | 4.682235831e+05 (×0.781767) |

and σ along the **`--b1-unpol`** axis at fixed `--unpol-sf` agrees to
**1–2 ulp**: over the five same-`unpol_sf` pairs the largest relative
difference is **2.5e−16**, and one pair (`toy` vs `mstw` at `--unpol-sf toy`)
is byte-equal. The tensor-thirds average cancels the b₁ term to **rounding,
not identically**, so this is stated as "to 2 ulp" and never as "bit for
bit". (`InclusiveSampler::cell_xsec_pb`, the SPIN-BLIND cell cross section
that `B1UnpolSource`'s own header calls bit-identical under this flag, is a
*different* quantity from the plan-weighted `sigma_pb` measured here.)

**A_zz moves on BOTH axes**, because it is the ratio of the two — at
x = 0.10, Q² = 2.5, y = 0.5, θ_m = 0:

| `b1_unpol` \ `unpol_sf` | `toy` | `ct18nlo` | `mstw` |
|---|---|---|---|
| `toy` | +1.100752e−05 | — | — |
| `ct18nlo` | +2.445101e−05 | +2.191100e−05 | +2.458518e−05 |
| `mstw` | +2.033932e−05 | +1.822644e−05 | +2.045092e−05 |

Over the seven legal cells A_zz at that point spans **+1.100752e−05 …
+2.458518e−05, a factor 2.23** — and that is *before* the mandatory ±100 %
band, which takes the same quantity to `[0, 2×]` of each entry.

---

## A4.3 The four terms — the band, x·b₁ per nucleon, ⁶Li, Q² = 2.5 GeV²

Default options throughout (raw digitized CDKS theory-1 b₁ᵈ, `r_sigma_lt`,
κ = 1, `norm_target` = `VMC_N_ALPHA_D_LI6`, w_CG = 1/10, w_αd = 1, the
2400/2400/3200 y grid). `LI6_B1_PER_NUCLEON` = 2/6 on the three deuteron
terms and 4/6 on the struck-α term.

**Terms (1) and (3) are bit-identical across all three backends** — verified as
doubles at all five x — because they convolve **b₁ᵈ**, not F₁. Only the
orbital sector moves. They are printed once.

| x | (1) embedded d, S | (3) CG D-wave |
|---|---|---|
| 0.05 | +1.568562e−6 | +3.097931e−9 |
| 0.10 | −4.692632e−6 | −9.296474e−9 |
| 0.20 | −2.366420e−5 | −4.673602e−8 |
| 0.30 | −4.496665e−5 | −8.834207e−8 |
| 0.50 | +4.244833e−5 | +8.564169e−8 |

### The orbital terms, which are the whole of the band

| x | | (2d) struck d | (2α) struck α | **total x·b₁** |
|---|---|---|---|---|
| 0.05 | `toy` | +1.383779e−6 | +6.952760e−7 | **+3.650716e−6** |
| | `ct18nlo` | −1.675024e−8 | +1.984903e−8 | **+1.574759e−6** |
| | `mstw` | +2.839977e−7 | +1.425989e−7 | **+1.998257e−6** |
| 0.10 | `toy` | +1.267261e−6 | +6.380547e−7 | **−2.796612e−6** |
| | `ct18nlo` | −1.004325e−6 | −5.058663e−7 | **−6.212120e−6** |
| | `mstw` | −3.351524e−7 | −1.304059e−7 | **−5.167487e−6** |
| 0.20 | `toy` | +1.666861e−6 | +8.388558e−7 | **−2.120522e−5** |
| | `ct18nlo` | −2.878556e−6 | −1.442155e−6 | **−2.803165e−5** |
| | `mstw` | −2.347871e−6 | −1.120821e−6 | **−2.717963e−5** |
| 0.30 | `toy` | +2.213341e−6 | +1.109171e−6 | **−4.173248e−5** |
| | `ct18nlo` | −4.394021e−6 | −2.250548e−6 | **−5.169956e−5** |
| | `mstw` | −5.432455e−6 | −2.761567e−6 | **−5.324901e−5** |
| 0.50 | `toy` | +6.163151e−6 | +3.042217e−6 | **+5.173934e−5** |
| | `ct18nlo` | +7.670734e−6 | +3.907924e−6 | **+5.411263e−5** |
| | `mstw` | −1.878181e−8 | −2.456684e−7 | **+4.226952e−5** |

The four terms sum to `b1` to **≤ 3.1e−16 relative** on every one of the
fifteen rows (the algebraic identity is configuration-independent, as
`phase_D_numbers.md` said it was).

### THE QUOTABLE BAND — the envelope, and then the mandatory ±100 % on top

| x | `--b1-unpol` envelope of x·b₁ | **× {0, 1, 2} — what may actually be quoted** |
|---|---|---|
| 0.05 | [+1.574759e−6, +3.650716e−6] | **[0, +7.301431e−6]** |
| 0.10 | [−6.212120e−6, −2.796612e−6] | **[−1.242424e−5, 0]** |
| 0.20 | [−2.803165e−5, −2.120522e−5] | **[−5.606329e−5, 0]** |
| 0.30 | [−5.324901e−5, −4.173248e−5] | **[−1.064980e−4, 0]** |
| 0.50 | [+4.226952e−5, +5.411263e−5] | **[0, +1.082253e−4]** |

The right-hand column is the honest one. The unpolarised-backend spread
(max/min over the three backends) is **×2.318 / ×2.221 / ×1.322 / ×1.276 /
×1.280** at x = 0.05 / 0.10 / 0.20 / 0.30 / 0.50, while the ±100 % band is a
factor ∞ downward (it reaches zero) and exactly ×2 upward.

**So at the two lowest x the PDF systematic is NOT contained inside the A > 2
band, and at the other three it is.** Taking the narrowest backend's own
`{0, 1, 2}` band, the widest backend's central value falls **outside** it at
x = 0.05 (×2.318 > 2) and at x = 0.10 (×2.221 > 2), and inside it at x = 0.20,
0.30 and 0.50. This is measured, and it is the reason both bands have to be
quoted rather than the wider one alone: at x ≤ 0.1 running `{0, 1, 2}` on one
backend does **not** bracket the other backends.

### The unpolarised-backend ratios, for the record

`b1(x, 2.5)` relative to the `toy` default:

| | x = 0.05 | 0.10 | 0.20 | 0.30 | 0.50 |
|---|---|---|---|---|---|
| `ct18nlo`/`toy` | 0.431356 | 2.221302 | 1.321922 | 1.238833 | 1.045870 |
| `mstw`/`toy` | 0.547360 | **1.847766** | 1.281742 | **1.275961** | **0.816971** |

The three bold entries are the **1.848 / 1.276 / 0.817** that
`OPEN_ITEMS_SOLUTIONS.md` §10, `USAGE.md` §2a, `PHYSICS_CHANNELS.md`,
`python/lipolgen/cli.py` and `python/bindings.cpp` all quote; they reproduce
here exactly. Note the two ends this run adds: at **x = 0.05 both real fits
are a factor ~2 BELOW the toy** (0.43 and 0.55), so the spread is not
one-sided, and it is **not monotone in x on either backend**.

### The reading that does NOT survive the band

`phase_D_numbers.md` and §10 both read the toy column as:

> *"Terms (2d)+(2α) are 7–133 % of term (1) and have the **opposite sign** over
> the whole window … They very nearly cancel term (1) at x ≈ 0.1 and dominate
> below. This is the ⁶Li quadrupole puzzle showing up in b₁."*

**The signed ratio [(2d)+(2α)]/(1), measured on all three backends:**

| backend | x = 0.05 | 0.10 | 0.20 | 0.30 | 0.50 | 0.70 |
|---|---|---|---|---|---|---|
| `toy` | +1.325453 | **−0.406023** | **−0.105886** | **−0.073888** | +0.216861 | +0.092585 |
| `ct18nlo` | +0.001976 | +0.321822 | +0.182584 | +0.147767 | +0.272771 | +0.166692 |
| `mstw` | +0.271967 | +0.099210 | +0.146580 | +0.182224 | −0.006230 | +0.165628 |

**The opposite sign is a `toy` feature at three of the five table x, and it
survives on neither real global fit.** On `ct18nlo` the orbital sector has the
**same** sign as term (1) at every x in the window; on `mstw` at four of five
(the −0.6 % at x = 0.50 is a near-zero, not an opposition). So the sentence
"they very nearly cancel term (1) at x ≈ 0.1" is true **only** on the toy
backend: on `ct18nlo` and `mstw` at x = 0.10 the orbital sector **adds to**
term (1), so |x·b₁| **exceeds** |term (1)| there (×1.324 on `ct18nlo`, ×1.101
on `mstw`) where on `toy` it is **0.596** of it. The totals are ×2.221 and
×1.848 the toy's — the latter being the same **1.848** already published at
that x, now with its mechanism named.
**What survives** is the statement the band is actually for: the orbital sector
is 0.2 %–133 % of term (1) and is the term the whole unpolarised-backend
spread lives in, so it is the quantitative reason the ±100 % band is mandatory.
**Its SIGN relative to term (1) may not be quoted without the backend.**

### (2α)/(2d) — the analytic 0.5064, and where it stops being readable

| backend | x = 0.05 | 0.10 | 0.20 | 0.30 | 0.50 | 0.70 |
|---|---|---|---|---|---|---|
| `toy` | 0.502447 | 0.503491 | 0.503255 | 0.501130 | 0.493614 | 0.502937 |
| `ct18nlo` | **−1.184999** | 0.503688 | 0.500999 | 0.512184 | 0.509459 | 0.508424 |
| `mstw` | 0.502113 | **0.389094** | 0.477378 | 0.508346 | **13.080123** | 0.502902 |

Against the analytic **2(M_d/M_α)² = 0.5064**. The three bold entries are
**not** a broken identity and must not be quoted as a spread: they are ratios
of two quantities that both pass through zero there — (2d) is −1.675024e−8 at
`ct18nlo`/x = 0.05 and −1.878181e−8 at `mstw`/x = 0.50, i.e. **×0.0121 and
×0.0030** of their `toy` values (+1.383779e−6 and +6.163151e−6). Away from
those crossings the ratio holds 0.477–0.512 on all
three backends, i.e. the T15 identity is configuration-robust wherever it is
readable.

---

## A4.4 Every other `phase_D_numbers.md` table, banded

### SD / DD split of the orbital sector (x·b₁)

| x | backend | SD | DD | SD/DD |
|---|---|---|---|---|
| 0.10 | `toy` | +1.753696e−6 | +1.516195e−7 | 11.5664 |
| | `ct18nlo` | −1.391140e−6 | −1.190509e−7 | 11.6853 |
| | `mstw` | −4.289026e−7 | −3.665573e−8 | 11.7008 |
| 0.30 | `toy` | +3.057410e−6 | +2.651024e−7 | 11.5329 |
| | `ct18nlo` | −6.117354e−6 | −5.272146e−7 | 11.6032 |
| | `mstw` | −7.541047e−6 | −6.529753e−7 | 11.5487 |
| 0.50 | `toy` | +8.468950e−6 | +7.364175e−7 | 11.5002 |
| | `ct18nlo` | +1.063033e−5 | +9.483307e−7 | 11.2095 |
| | `mstw` | −2.503738e−7 | −1.407633e−8 | **17.7869** |

SD dominates DD by **11.21–11.70×** on every readable entry — a range **0.49
wide** against the published "11.50–11.57×", which was the `toy` column's own
**0.07**. The 17.79 at `mstw`/x = 0.50 is the same near-zero crossing as
above, not a different physics.

### Close–Kumano integrals — REPORTED, not enforced

| model | ∫b₁ dx, [0.01, 1.2], 241 points |
|---|---|
| `Li6ConvolutionB1`, `--b1-unpol toy` | +1.458804e−4 |
| `Li6ConvolutionB1`, `--b1-unpol ct18nlo` | +1.394021e−4 |
| `Li6ConvolutionB1`, `--b1-unpol mstw` | +1.387049e−4 |
| **band** | **[+1.387049e−4, +1.458804e−4]** — a **5.2 %** spread |
| `Li6B1(MillerB1)`, same grid | +9.124294e−4 |
| `Li6B1(CdksB1)`, same grid | +1.385959e−4 |
| `close_kumano_integral(true)` (digitized CDKS, its own grid) | +4.592003e−4 |
| `close_kumano_integral(false)` (Miller, its own grid) | +5.914740e−3 |

The integral is the **tightest** quantity in this document across the axis
(**5.17 %**, against ×1.28–×2.32 pointwise) — the orbital sector's sign change
in x partly cancels under ∫dx. None is zero; none is enforced.

### The mandatory band rows — `--b1-band-scale 0 / 1 / 2` (x·b₁ at x = 0.30)

| backend | scale 0 | scale 1 | scale 2 | scale 2 == 2 × scale 1? |
|---|---|---|---|---|
| `toy` | 0 (exactly) | −4.173248e−5 | −8.346496e−5 | **yes, as doubles** |
| `ct18nlo` | 0 (exactly) | −5.169956e−5 | −1.033991e−4 | **yes, as doubles** |
| `mstw` | 0 (exactly) | −5.324901e−5 | −1.064980e−4 | **yes, as doubles** |

The band's exact linearity is algebraic and survives the regeneration, as
predicted.

### Term-(2) knob rows — `--b1-alpha-d-dwave-weight` 0 / 1 / 2

x·b₁ at x = 0.10 / 0.30 / 0.50:

| backend | w = 0 | w = 1 | w = 2 |
|---|---|---|---|
| `toy` | −4.701928e−6 / −4.505499e−5 / +4.253397e−5 | −2.796612e−6 / −4.173248e−5 / +5.173934e−5 | −8.912966e−7 / −3.840997e−5 / +6.094470e−5 |
| `ct18nlo` | *(identical)* | −6.212120e−6 / −5.169956e−5 / +5.411263e−5 | −7.722311e−6 / −5.834413e−5 / +6.569129e−5 |
| `mstw` | *(identical)* | −5.167487e−6 / −5.324901e−5 / +4.226952e−5 | −5.633045e−6 / −6.144303e−5 / +4.200507e−5 |

**The w = 0 row is bit-identical on all three backends** (verified as doubles)
— which is the same fact as A4.2's, seen from the other side: w = 0 removes the
orbital terms and with them the entire `--b1-unpol` dependence. The knob stays
exactly linear on every backend. The published reading *"at x = 0.10 the
orbital term is 70 % of the total's magnitude, so w = 2 nearly zeroes b₁
there"* is **toy-only**: on `ct18nlo` and `mstw`, w = 2 moves x·b₁(0.10) by
+24 % and +9 % **away** from zero.

### `finite_q_delta` off (the ⁶Li default) / on

| x | backend | κ = 1 (default) | κ = √(1+γ²) | ratio κ/κ=1 |
|---|---|---|---|---|
| 0.10 | `toy` | −2.796612e−6 | −2.769979e−6 | 0.990477 |
| | `ct18nlo` | −6.212120e−6 | −6.233773e−6 | 1.003486 |
| | `mstw` | −5.167487e−6 | −5.175688e−6 | 1.001587 |
| 0.30 | `toy` | −4.173248e−5 | −4.125623e−5 | 0.988588 |
| | `ct18nlo` | −5.169956e−5 | −5.248069e−5 | 1.015109 |
| | `mstw` | −5.324901e−5 | −5.423648e−5 | 1.018544 |
| 0.50 | `toy` | +5.173934e−5 | +5.553091e−5 | 1.073282 |
| | `ct18nlo` | +5.411263e−5 | +5.875855e−5 | 1.085857 |
| | `mstw` | +4.226952e−5 | +4.293734e−5 | 1.015799 |

**Across the whole band the switch is −1.1 % to +8.6 %**, against 1.57–1.69 on
the A = 2 gate. The published "far smaller than on the A = 2 gate" reading is
**unchanged and is now a band statement**: this is why the two objects default
differently (design A10) — one cached f(y) table against a ~0.24 s rebuild per
x, for at most 8.6 %.

### Densities — MEASURED to be independent of the axis

| quantity | `toy` | `ct18nlo` | `mstw` |
|---|---|---|---|
| `norm()` (renormalised), struck d / struck α | 0.8194810000 / 0.8194810000 | *identical* | *identical* |
| `mean_y()` = ⟨z²⟩/⟨z⟩, struck d / struck α | 0.9997871376 / 0.9984638289 | *identical* | *identical* |
| `p_d()` / `p_d_momentum()` | 0.0193299187 / 0.0193516197 | *identical* | *identical* |
| `y_max`, struck d / struck α | 1.9928570 / 1.2512039 | *identical* | *identical* |

All eight quantities agree to every printed digit on all three backends, as
they must — `LightConeDensities` is built from the VMC α–d tables and never
sees an `UnpolSF`. `VMC_P_D_LI6` = 0.019354933,
`VMC_N_ALPHA_D_LI6` = 0.8194810. **These rows of `phase_D_numbers.md` were
never a ToyF2 measurement and did not need regenerating**; that is recorded
here so nobody re-runs them looking for a move.

### T14 — the angular-average truncation

| x | P₂ remainder / (1) | P₄ remainder / (1) |
|---|---|---|
| 0.05 | −2.312802e−5 | +5.402826e−7 |
| 0.10 | −8.730694e−6 | −9.662514e−6 |
| 0.20 | −1.702969e−6 | −1.048166e−6 |
| 0.30 | +1.410904e−5 | −4.718219e−7 |
| 0.50 | −5.105656e−5 | −1.360977e−5 |

**Bit-identical on all three backends** (verified as doubles): both remainders
and their denominator, term (1), convolve b₁ᵈ rather than F₁, so the ratio
carries no unpolarised backend at all. The "below 6e−5 of term (1) everywhere
in the window" reading is therefore configuration-independent, not a toy
measurement — again recorded so it is not re-run.

### The ±5 % N_αd systematic and the Q4 knob (x·b₁ at x = 0.30)

| backend | 0.95 × N_αd | N_αd (default) | 1.05 × N_αd | Q4 (`use_spectroscopic_factor = false`) |
|---|---|---|---|---|
| `toy` | −3.964585e−5 | −4.173248e−5 | −4.381910e−5 | −5.092550e−5 |
| `ct18nlo` | −4.911458e−5 | −5.169956e−5 | −5.428454e−5 | −6.308817e−5 |
| `mstw` | −5.058656e−5 | −5.324901e−5 | −5.591146e−5 | −6.497895e−5 |

b₁ is **exactly linear** in `norm_target` on every backend (0.950000 /
1.050000 as doubles), and Q4 is **×1.220285 = 1/N_αd** on every backend, to
six digits. Both systematics are therefore **multiplicative and
configuration-independent**: ±5 % flat, and ×1.2203 for Q4, whatever the
unpolarised backend. This is the one place where a published ⁶Li number
genuinely did not need the rerun, and it is because the quantity is a ratio.

### Alternative deuteron camp — `deuteron_b1 = MillerB1` (x·b₁)

| backend | x = 0.10 | x = 0.30 | x = 0.50 |
|---|---|---|---|
| `toy` | +3.610434e−4 | **+6.230330e−6** | −1.495634e−4 |
| `ct18nlo` | +3.576279e−4 | **−3.736752e−6** | −1.471901e−4 |
| `mstw` | +3.586725e−4 | **−5.286205e−6** | −1.590332e−4 |

At x = 0.10 and 0.50 the camp swamps the backend choice (0.9 % and 8.0 %
spread). **At x = 0.30 the band spans zero and the sign is a band edge**:
+6.23e−6 on `toy` against −3.74e−6 and −5.29e−6 on the two real fits. The
published +6.2303e−6 stands only as a band edge, and this is a legend row for
a plot in any case — Miller and CDKS are different *camps* for b₁ᵈ and the
default does not change.

---

## A4.5 What this changes in the standing claims

| claim, as published | status after the band |
|---|---|
| the four terms sum to `b1` to 1e−12 | **holds** on all three backends (≤ 3.1e−16) |
| the band's exact linearity (scale 2 == 2 × scale 1) | **holds** bit for bit on all three |
| term (3) is 0.197 % of term (1) at every x | **holds and is now known to be exact** across the axis — both terms are backend-independent |
| (2α)/(2d) ≈ 0.5064 | **holds**, 0.477–0.512 wherever readable; three entries are near-zero crossings and are not a spread |
| SD/DD ≈ 11.5 | **widens** to 11.21–11.70 |
| ±5 % N_αd, ×1.2203 Q4 | **hold exactly** on all three — multiplicative |
| T14 remainders ≤ 6e−5 of term (1) | **holds and is backend-independent** (bit-identical) |
| the densities | **backend-independent** (identical) |
| `finite_q_delta` is −1 % to +7 % | **restated as −1.1 % to +8.6 %** over the band |
| ∫b₁ dx = +1.4588e−4 | **becomes a band**, [+1.387049e−4, +1.458804e−4] |
| every x·b₁ **value** | **becomes a band**; see A4.3 |
| **"terms (2d)+(2α) have the opposite sign to term (1)"** | **DOES NOT SURVIVE.** True on `toy` at x = 0.10/0.20/0.30 only; **false at every x on `ct18nlo`** and at four of five on `mstw` |
| **"w = 2 nearly zeroes b₁ at x = 0.10"** | **DOES NOT SURVIVE.** Toy-only; on the two real fits w = 2 moves b₁(0.10) *away* from zero by +24 % / +9 % |

---

## A4.6 What is pinned, and where

`python/tests/test_li6_unpol_band.py` pins **one row of the band** — the
x = 0.30, Q² = 2.5 row of A4.3, all three backends, at rtol 1e-12 — plus the
two structural facts the band rests on (terms (1)/(3) bit-identical across the
axis; the `--unpol-sf` axis exactly flat on x·b₁ through the pipeline, and the
two refused cells refused by name). The `mstw` and `ct18nlo` rows skip loudly
in a build without the matching optional tier, exactly as the gate's own rows
do.

## A4.7 How to reproduce

```bash
cd /path/to/LiPolGen && source env.sh
cmake --build build -j
python -m pytest python/tests/test_li6_unpol_band.py -q
```

and the two snippets of §A4.1 for the full tables. **Nothing in this document
is written by a script into `docs/`**; `validation/b1_li6_table.py --write`
still emits only its own `run_2026-09-02/phase_D_regenerated.md` on the
default backend and is untouched by this section.

## A4.8 The other half of A4 — registry row 6 / §B11, republished BESIDE the original

Task A4 carries a second, unrelated item: `run_2026-09-02/phase_C_numbers.md`
§8.1c/§8.2's `A_zz(Born)` column is **×3.253983 too large**, and the registry's
options were "leave with the note" or "republish". **The thing that takes no
decision away was done**: the originals stay, labelled with what computed
them, and the recomputation is a dated table next to each —
`phase_C_numbers.md` **§8.1c-corr** (all eight (x, Q²) points × both C0
edges), **§8.2-corr** (the five-row band with τ, `w_hi − 1` and the
half-width) and **§8.3-corr** (the six run-level rows that inherit the
column); `OPEN_ITEMS_SOLUTIONS.md` §9 carries the same pairs at the same date.
Both registry options remain open and (ii) now costs one deletion.

**The factor is exact and it is named.** `published ÷ shipped = 3.253983147` at
every point, which is `1 / (LI6_B1_RANK2_TRANSFER × LI6_B1_PER_NUCLEON)` =
`1/(0.921947 × 2/6)` — the two factors `Li6B1` applies to the deuteron table
and `toy_b1` does not. It is exact because the shipped `b2` is the
Callan–Gross `2·x·b₁` (measured `tables().b2 == 2*x*tables().b1` as doubles),
which makes `azz` linear in b₁. **`ΔA_zz` does not divide by it**, because
`ΔA_zz = [2 r_T − A_zz r_U]/(1 + r_U)` has an `A_zz`-independent term.

**Four findings, each measured** (the numbers and the recipe live in
`phase_C_numbers.md` §8.1c-corr, not here, so no physics number is defined
twice):

1. `ΔA_zz` changes **sign at x = 0.01 on both C0 edges at all three Q²** —
   six entries, where the 2026-09-04 note recorded one.
2. Which knob **spreads** move follows one rule: a knob that leaves `r_U`
   alone has an `A_zz`-independent spread. `fq_scale` and `tail_tensor_scale`
   keep theirs; `qe_suppression` and `qe_kf_gev` shrink ×0.31; `c0_shape`
   **widens** +15 %.
3. The RC budget's **order changes at second place** once §9's already-
   corrected `qe_tensor_scale` row is compared on the same footing as the four
   that were not.
4. §9's conclusions **1 and 2 do not depend on the column** (re-measured
   identical to every printed digit); **3 does and stands** — ×767.6, the 768
   already on record.

One reproduction exception is recorded rather than smoothed: three of the
eight `vmc-ft` `ΔA_zz` entries re-measure 2.4–6.3 × 10⁻⁵ relative away from
their published values (fifth significant figure), deterministically in this
tree, cause not established; every `ho` row — which carries all the
conclusions — is exact.

---

## A4.9 End state, re-measured after everything above

Nothing in this section changed a source file; both halves of A4 are
documentation plus one new pytest file. Re-measured on the finished tree:

* **C++ (doctest): 404 cases / 17 240 384 assertions / 1 skipped / 0 failed** —
  unchanged from the section's own baseline.
* **pytest: 1004 passed / 127 skipped** — the baseline's 992 plus the **12**
  of `python/tests/test_li6_unpol_band.py`; the 127 skips are the pre-existing
  refused-base matrix cells and are unchanged.
* **Docs gate: 1220 references (strict) / 95 ranges / 7 external / 0 broken /
  6 allow-listed**, +19/115 (`SPIN32_FINITE_GAMMA.md`), +8 external
  (`PYTHIA_BRIDGE.md`). The `USAGE.md` §2a rewrite pushed six cited **ranges**
  down the file; they were re-anchored with `--fix` (block moved, content
  untouched — `docs/USAGE.md:1810-2001 → :1860-2051`, `:2029-2204 → :2079-2254`,
  `:1352-1358 → :1402-1408`, `:1378-1386 → :1428-1436`, the last two cited
  twice each). **No `--record-ranges` was needed or used**: no cited block was
  edited.
* **SPDX: 101/101.**
* **`git diff --stat validation/reference/` is empty** — no rtol-1e−12
  reference moved, and none could: this section writes documentation and reads
  the library.
