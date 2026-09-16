# Run 2026-09-06 — Phase A measurements

Baseline this phase started from, re-measured before anything moved (commit
`91e48b9`, tree otherwise clean): **401 doctest cases / 17 240 286 assertions
/ 1 skipped / 0 failed**; **926 pytest passed / 112 skipped**; docs gate
**1206 references (strict) / 96 ranges / 7 external / 0 broken / 6
allow-listed**, +19/115 (`SPIN32_FINITE_GAMMA.md`), +8 external
(`PYTHIA_BRIDGE.md`); SPDX **101/101**.

---

## A1 — registry row 3, option (iii): ONE R hook in BOTH halves of the ⁶Li tensor weight

**What this section prices, and what it does not decide.** `STATUS.md` decision
row 3 (`../run_2026-09-03/AUTHOR_DECISIONS.md` §B7) asks which R = σ_L/σ_T the
opt-in ⁶Li convolution backend should use, and offers three answers: **(i)**
`r_sigma_lt`, the status quo; **(ii)** `r1998`, the R the A = 2 gate was
validated with; **(iii)** *one* R hook threaded through
`default_inclusive_kernel` into both the convolution and the kernel, which the
run itself called "the better third option, not taken here" and left
**unpriced**. (iii) is now implemented, **opt-in, default unchanged**, and
priced here. **Nothing is decided: no default moved, and every reference JSON
is untouched at rtol 1e-12.**

### A1.1 Why a different R in the numerator does not cancel

The shipped observable is a **ratio**. `azz` (`include/lipolgen/asymmetries.hpp`)
is

    A_zz = TENSOR_LL_SIGN · (2/3) · K / D_φ ,   TENSOR_LL_SIGN = −1,

with **K = b₁ + (1−y)/(x y²) b₂** and **D_φ = F₁ + (1−y)/(x y²) F₂**
(`InclusiveKernel::tensor_kernel` and `::dphi`). At θ_m = 0 the geometry factor
is 1, so **A_zz is exactly −(2/3) of the weight K/D_φ** — one number, not two,
and every percentage below is the same for both.

The numerator's R lives on `Li6ConvolutionOptions::r_func` and the
denominator's on `InclusiveKernel::Options::r_func`. Both are null by default,
and `resolve_r` (`include/lipolgen/sf.hpp`) turns a null hook into
`r_sigma_lt` — so the two halves agree today **by the coincidence of two
independent defaults**, and before this change the only way to move one was to
move it alone. That is option (ii), and the factor it fails to cancel is
(1 + r1998)/(1 + r_sigma_lt), tabulated in A1.3.

### A1.2 The wiring, and why the default cannot move

`PipelineConfig::r_source` (`RSource`; CLI `--r-source`, Python
`lipolgen.R_SOURCE`) has three values:

| value | what it installs | status in `knob_provenance` |
|---|---|---|
| `unset` (**default**) | **nothing** — the branch is not taken, both hooks stay null | `read` (another value gives another file) |
| `sigma-lt` | ONE `r_sigma_lt` object in `Li6ConvolutionOptions::r_func` **and** `InclusiveKernel::Options::r_func` | `not-read` — measured bit-identical to `unset` |
| `r1998` | ONE `r1998` object in both | `read` |

`unset` does not enter the wiring branch at all. That is deliberate: it makes
the default bit for bit **by construction**, not by an argument that an
explicit `r_sigma_lt` closure resolves to the same double as `resolve_r`'s null
branch. (It does — §A1.4 — but the default does not rest on it.)

**Reach.** The hook is installed by the `Li6Convolution` branch of
`default_inclusive_kernel` and nowhere else, so `PipelineConfig::validate()`
refuses anything but `unset` off `--b1-model li6-convolution` — which is
already inclusive-only, ⁶Li-only, and refused beside a caller-supplied
`kernel`. That one condition is the whole reach rule. `meta["r_source"]` is
written unconditionally, which is why a value that did not reach the rate is
refused rather than recorded. On `miller` and `cdks` the numerator has **no F₁
of its own** — Miller is a ratio model, CDKS a digitized column — so "one
shared R" is not a statement those branches can make.

**The asymmetry with `--b1-unpol`, and it is the point of (iii).**
`--b1-unpol` moves the numerator alone and leaves the spin-blind cell cross
section bit for bit. `--r-source` moves **both**: the kernel's `r_func` is
threaded to all four places the kernel needs R — F₁ in `NuclearF2::f1a`, F_L in
`dsigma_dx_dq2`, D(y) in `depolarization_d`, and the default `ToyG1`'s F₁ — so
the unpolarised rate moves with it (§A1.6).

### A1.3 R itself, and the numerator (x·b₁)

⁶Li, `--b1-model li6-convolution`, shipped defaults (raw digitized CDKS b₁ᵈ,
`ToyF2`, κ = 1, N_αd = 0.819481, band 1). x·b₁ is **y-independent**.

| x | Q² | `r_sigma_lt` | `r1998` | (1+r1998)/(1+r_sigma_lt) | x·b₁ under (i) and (iii) `sigma-lt` | x·b₁ under (ii) and (iii) `r1998` | ratio |
|---|---|---|---|---|---|---|---|
| 0.05 | 2.5 | 0.171429 | 0.300539 | 1.110217 | +3.650716e−6 | +3.449387e−6 | 0.944852 |
| 0.05 | 5.0 | 0.163636 | 0.212035 | 1.041593 | +4.875303e−6 | +4.724903e−6 | 0.969151 |
| 0.10 | 2.5 | 0.171429 | 0.273940 | 1.087510 | -2.796612e−6 | -3.484857e−6 | 1.246099 |
| 0.10 | 5.0 | 0.163636 | 0.184434 | 1.017873 | -1.777795e−6 | -2.570765e−6 | 1.446042 |
| 0.30 | 2.5 | 0.171429 | 0.217570 | 1.039389 | -4.173248e−5 | -4.284379e−5 | 1.026629 |
| 0.30 | 5.0 | 0.163636 | 0.130957 | 0.971916 | -4.031080e−5 | -4.120370e−5 | 1.022151 |

**(ii) and (iii) `r1998` give bit-identical b₁ and b₂** — both hand `r1998` to
`Li6ConvolutionOptions`, so the numerator is the same object's output and every
difference in A1.4 is the **denominator's**. The Q² = 2.5 column reproduces
§B7(d)'s −5.5 % / +24.6 % / +2.7 % at x = 0.05 / 0.10 / 0.30 exactly. Note the
**sign flip at x = 0.30, Q² = 5**: there `r1998` = 0.1310 is *below*
`r_sigma_lt` = 0.1636 and the mismatch factor is 0.9719, not > 1 — the
registry's "1.088 / 1.039" is a Q² = 2.5 statement and does not carry to
Q² = 5.

### A1.4 The tensor observables at the standard points

x = 0.05 / 0.10 / 0.30 × Q² = 2.5 / 5 GeV², at **y = 0.5** — this tree's
convention for structure-function-level asymmetry tables
(`../run_2026-09-03/phase_D_numbers.md` §, `phase_D_sf_injection.md`). y is a
free variable here and is **not** the beam kinematics' y; §A1.5 is the window.
A_zz = −(2/3) K/D_φ, so its percentage column is the weight's and is not
repeated.

| x | Q² | K/D_φ (i) | K/D_φ (ii) | K/D_φ (iii) `r1998` | (ii)−(i) | (iii)−(i) | A_zz (i) | A_zz (iii) | cos 2φ (i) | cos 2φ (iii) | (iii)−(i) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.05 | 2.5 | +1.837372e−5 | +1.736046e−5 | +1.766896e−5 | -5.5148 % | -3.8357 % | -1.224915e−5 | -1.177931e−5 | -1.166354e−3 | -1.069234e−3 | -8.3268 % |
| 0.05 | 5.0 | +2.318077e−5 | +2.246565e−5 | +2.262543e−5 | -3.0849 % | -2.3957 % | -1.545385e−5 | -1.508362e−5 | -1.172783e−3 | -1.133960e−3 | -3.3104 % |
| 0.10 | 2.5 | -1.651128e−5 | -2.057469e−5 | -2.087005e−5 | +24.6099 % | +26.3988 % | +1.100752e−5 | +1.391337e−5 | -1.156685e−3 | -1.078877e−3 | -6.7268 % |
| 0.10 | 5.0 | -1.017691e−5 | -1.471623e−5 | -1.476207e−5 | +44.6042 % | +45.0546 % | +6.784604e−6 | +9.841380e−6 | -1.163061e−3 | -1.146198e−3 | -1.4499 % |
| 0.30 | 2.5 | -4.695082e−4 | -4.820109e−4 | -4.852452e−4 | +2.6629 % | +3.3518 % | +3.130055e−4 | +3.234968e−4 | -5.885366e−4 | -5.700326e−4 | -3.1441 % |
| 0.30 | 5.0 | -4.512556e−4 | -4.612512e−4 | -4.589061e−4 | +2.2151 % | +1.6954 % | +3.008371e−4 | +3.059374e−4 | -5.917807e−4 | -6.057850e−4 | +2.3665 % |

**(iii) with `sigma-lt` is bit-identical to (i)** — 162 values (b₁, b₂, F₁,
F₂, Δ, D_φ, K/D_φ, A_zz, cos 2φ at 3 x × 2 Q² × 3 y), every one equal as a
double, so it does not get a column. That is the wiring's own test: the shared
hook reaches both halves and changes nothing when it names the value they
already had.

**Read this table as three statements.**

1. **The cos 2φ amplitude separates the two options cleanly.** Its numerator is
   Δ = `toy_delta_gluon`, which carries **no R at all**, so option (ii) leaves
   it *exactly* unchanged (bit for bit, at every point) and (iii) moves it
   **purely through D_φ**: −8.33 % / −6.73 % / −3.14 % at Q² = 2.5 and
   −3.31 % / −1.45 % / +2.37 % at Q² = 5. Anyone extracting Δ from the cos 2φ
   amplitude is unaffected by (ii) and is not by (iii).
2. **On the weight the two options differ by the denominator's move alone**:
   at Q² = 2.5, (ii) gives −5.5148 % / +24.6099 % / +2.6629 % and (iii)
   −3.8357 % / +26.3988 % / +3.3518 %. The gap is **1.7 / 1.8 / 0.7 points** —
   the same order as the shift itself at x = 0.05 and x = 0.30, and small
   beside it at x = 0.10, where the +24.6 % is the near-cancellation of term
   (1) against the orbital terms and not a (1 + R) prefactor.
3. **(iii) does not make the shift smaller.** The registry's framing — the
   numerator's (1 + R) "no longer cancels the one in the denominator" — is
   right about the mechanism and would be wrong if read as "sharing the hook
   restores the cancellation". It restores **part** of it, and at x = 0.10 it
   makes the shift *larger*. The cancellation is only complete where D_φ is
   pure F₁, i.e. as y → 1 (§A1.5), and the ⁶Li b₁ terms are not proportional
   to F₁ anyway: they convolve F₁ᵈ against a density that integrates to zero
   and so respond to the **slope** of R.

### A1.5 The y window — the signature of (iii)

Under (ii) the weight's shift is **y-independent**. D_φ does not move at all,
and at the default b₂ = 2x·b₁ the numerator is K = b₁·[1 + 2(1−y)/y²] — the
whole y factor is a multiplier on b₁ — so the (ii)/(i) ratio of K/D_φ is
exactly b₁′/b₁ at every y. Measured: the three y agree to 2 ulp. Under (iii)
D_φ moves as well, and its F₁ and F₂ terms carry different powers of y (only
F₁ has an R in it), so the shift runs with y.

y = 0.1:

| x | Q² | K/D_φ (i) | K/D_φ (ii) | K/D_φ (iii) `r1998` | (ii)−(i) | (iii)−(i) | A_zz (i) | A_zz (iii) | cos 2φ (i) | cos 2φ (iii) | (iii)−(i) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.05 | 2.5 | +1.785039e−5 | +1.686598e−5 | +1.687389e−5 | -5.5148 % | -5.4705 % | -1.190026e−5 | -1.124926e−5 | -1.408591e−3 | -1.269348e−3 | -9.8853 % |
| 0.05 | 5.0 | +2.254633e−5 | +2.185078e−5 | +2.185493e−5 | -3.0849 % | -3.0666 % | -1.503088e−5 | -1.456995e−5 | -1.417978e−3 | -1.361614e−3 | -3.9750 % |
| 0.10 | 2.5 | -1.604099e−5 | -1.998866e−5 | -1.999626e−5 | +24.6099 % | +24.6573 % | +1.069399e−5 | +1.333084e−5 | -1.396913e−3 | -1.284994e−3 | -8.0119 % |
| 0.10 | 5.0 | -9.898371e−6 | -1.431346e−5 | -1.431465e−5 | +44.6042 % | +44.6162 % | +6.598914e−6 | +9.543101e−6 | -1.406223e−3 | -1.381647e−3 | -1.7477 % |
| 0.30 | 2.5 | -4.561353e−4 | -4.682819e−4 | -4.683657e−4 | +2.6629 % | +2.6813 % | +3.040902e−4 | +3.122438e−4 | -7.107681e−4 | -6.839549e−4 | -3.7724 % |
| 0.30 | 5.0 | -4.389051e−4 | -4.486271e−4 | -4.485655e−4 | +2.2151 % | +2.2010 % | +2.926034e−4 | +2.990437e−4 | -7.155051e−4 | -7.360790e−4 | +2.8754 % |

y = 0.9:

| x | Q² | K/D_φ (i) | K/D_φ (ii) | K/D_φ (iii) `r1998` | (ii)−(i) | (iii)−(i) | A_zz (i) | A_zz (iii) | cos 2φ (i) | cos 2φ (iii) | (iii)−(i) |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.05 | 2.5 | +2.020758e−5 | +1.909318e−5 | +2.068606e−5 | -5.5148 % | +2.3678 % | -1.347172e−5 | -1.379070e−5 | -3.175163e−4 | -3.098545e−4 | -2.4130 % |
| 0.05 | 5.0 | +2.539254e−5 | +2.460919e−5 | +2.539700e−5 | -3.0849 % | +0.0175 % | -1.692836e−5 | -1.693133e−5 | -3.179909e−4 | -3.150661e−4 | -0.9198 % |
| 0.10 | 2.5 | -1.815924e−5 | -2.262822e−5 | -2.413458e−5 | +24.6099 % | +32.9052 % | +1.210616e−5 | +1.608972e−5 | -3.148841e−4 | -3.088210e−4 | -1.9255 % |
| 0.10 | 5.0 | -1.114793e−5 | -1.612037e−5 | -1.634329e−5 | +44.6042 % | +46.6038 % | +7.431951e−6 | +1.089553e−5 | -3.153547e−4 | -3.141017e−4 | -0.3973 % |
| 0.30 | 2.5 | -5.163691e−4 | -5.301197e−4 | -5.461742e−4 | +2.6629 % | +5.7720 % | +3.442461e−4 | +3.641161e−4 | -1.602172e−4 | -1.588138e−4 | -0.8760 % |
| 0.30 | 5.0 | -4.943118e−4 | -5.052611e−4 | -4.941688e−4 | +2.2151 % | -0.0289 % | +3.295412e−4 | +3.294459e−4 | -1.604567e−4 | -1.614688e−4 | +0.6308 % |

**The window is wide and it straddles zero.** At x = 0.05, Q² = 2.5 the (iii)
shift runs **−5.4705 % at y = 0.1 → −3.8357 % at y = 0.5 → +2.3678 % at
y = 0.9**, against (ii)'s flat −5.5148 %; at x = 0.30, Q² = 5 it runs
+2.2010 % → +1.6954 % → −0.0289 %. **A single-y quotation of (iii) is not a
statement about a run**, which is the one thing (ii) does not suffer from. Any
published (iii) number must carry its y.

Where the y dependence comes from — F₁ and D_φ:

| x | Q² | F₁ (i) | F₁ (iii) | F₁ change | D_φ (i), y = 0.5 | D_φ (iii) | D_φ change, y = 0.1 / 0.5 / 0.9 |
|---|---|---|---|---|---|---|---|
| 0.05 | 2.5 | +3.494586e+0 | +3.147662e+0 | -9.9275 % | +1.986922e+1 | +1.952230e+1 | -0.0469 % / -1.7460 % / -7.7003 % |
| 0.05 | 5.0 | +3.719428e+0 | +3.570904e+0 | -3.9932 % | +2.103167e+1 | +2.088315e+1 | -0.0190 % / -0.7062 % / -3.1019 % |
| 0.10 | 2.5 | +1.489487e+0 | +1.369630e+0 | -8.0468 % | +8.468796e+0 | +8.348939e+0 | -0.0380 % / -1.4153 % / -6.2415 % |
| 0.10 | 5.0 | +1.544678e+0 | +1.517555e+0 | -1.7559 % | +8.734454e+0 | +8.707331e+0 | -0.0083 % / -0.3105 % / -1.3640 % |
| 0.30 | 2.5 | +2.605522e−1 | +2.506782e−1 | -3.7896 % | +1.481425e+0 | +1.471551e+0 | -0.0179 % / -0.6665 % / -2.9394 % |
| 0.30 | 5.0 | +2.632993e−1 | +2.709075e−1 | +2.8896 % | +1.488838e+0 | +1.496446e+0 | +0.0137 % / +0.5110 % / +2.2446 % |

At y → 0 the F₂ term of D_φ dominates and the denominator barely moves (≤ 0.05 %
here), so (iii) → (ii); at y → 1 the F₁ term dominates and the denominator
moves by the full (1 + R) ratio, which is where the cancellation the registry
describes actually happens.

### A1.6 Run level: (iii) moves the unpolarised rate, (ii) does not

⁶Li inclusive, config 1, `--b1-model li6-convolution`, `--x-max 0.95`,
`tensor-thirds` at P_z = 0.7, P_zz = 0.6, 2000 events, seed 7, one thread.

**The hash recipe, stated (2026-09-15).** As first published this table gave
two bare digests (`105e234556382ee6` / `2ace7fc1bc3ffa9e`) under the sentence
*"the hash is sha256 over every generated ndarray column"*, which does not fix
a digest: it names neither the column order nor whether the columns are hashed
separately or concatenated, and six plausible readings of it reproduce none of
the two. **They are replaced here by a recipe a reader can run** — one
`sha256`, columns taken in `sorted()` name order, each contributing
`name.encode("utf-8")` followed by `np.ascontiguousarray(col).tobytes()`,
truncated to the first 16 hex characters, over the 48 columns of the `--npz`
(`meta` excluded) — together with **the byte comparison the digests existed to
make**, which needs no recipe at all.

| `--r-source` | σ [pb] | σ(+1) ≈ σ(−1) [pb] (equal to 15 printed digits; 1 ulp apart at `unset`/`sigma-lt`, equal as doubles at `r1998`) | σ(0) [pb] | columns differing from `unset`, of 48 | sha256/16, recipe above | `meta["r_source"]` |
|---|---|---|---|---|---|---|
| `unset (default)` | 591846.161405217 | 591841.960163281 | 591854.563889088 | — | `da4e0c5cd2811390` | `unset` |
| `sigma-lt` | 591846.161405217 | 591841.960163281 | 591854.563889088 | **0 — byte-identical** | `da4e0c5cd2811390` | `not read at r_source = sigma-lt` |
| `r1998` | 587793.902124821 | 587789.915152379 | 587801.876069703 | **12** | `fc878f391d785aff` | `r1998` |

Every σ in this table reproduces to every printed digit on the committed tree
(`sigma_pb` 591846.1614052168 / 587793.9021248205; per-category
591841.960163281 / 591841.9601632811 / 591854.5638890883 and 587789.9151523793
/ 587789.9151523793 / 587801.8760697031); it is only the digests that were
unreproducible, and the conclusion they carried — `sigma-lt` byte-identical to
`unset`, `r1998` not — is re-measured directly above.

σ moves by **−0.684681 %** (ratio 0.993153188), and every per-category cross
section with it. Option (ii) **cannot** produce this: it never touches the
kernel's R, so a (ii) run has the same rate and only the tensor split moves.
That is the price of (iii) — the ratio's two halves agree, and the rate is no
longer the `r_sigma_lt` rate every published number was made on.

### A1.7 What (iii) costs, stated for the decision

* **Nothing at the default.** `unset` installs nothing; the branch is not
  taken. The default `miller` path, every reference JSON at rtol 1e-12, and
  the A = 2 gate (`DeuteronConvolutionB1`, which has its own `r_func` and is
  not touched) are all unmoved.
* **It removes the *reason* the row was decided the way it was, and does not
  remove the *effect*.** §B7's deciding argument is that "a different R in the
  numerator does not cancel". Under (iii) the two halves are the same choice
  by construction — but the shift against today is still −3.8 % to +26.4 % at
  y = 0.5, Q² = 2.5, and it is **y-dependent**, which (ii)'s is not.
* **It costs one new run-surface axis** with three values, one refusal rule,
  one `meta` key, one banner line, and a knob-provenance row whose middle
  value is `not-read`.
* **It buys the only configuration in which the tensor weight's numerator and
  denominator provably cannot disagree**, and it makes `--r-source r1998` the
  first way the run surface can put a *fitted* R anywhere: until now "R is
  hard-locked to `r_sigma_lt` and no CLI flag supplies one" was true of the
  whole program.
* **It does not settle (i) against (ii).** Choosing (iii) still requires
  naming an R. `sigma-lt` reproduces today; `r1998` is CDKS's and the gate's.
  The row therefore becomes: *(iii) with which R?* — and the numbers above are
  what that costs relative to each of (i) and (ii).

### A1.8 How to reproduce

```
cd LiPolGen && source env.sh && cmake --build build -j
build/lipolgen_tests -tc="*T18*r-source*"                            # 1 case, 58 assertions
python -m pytest python/tests/test_b1_model.py -k r_source -q        # 6 passed
python -m pytest python/tests/test_knob_provenance.py -k r_source -q # 9 passed, 12 skipped
```

The 12 skips are the ⁷Li, deuteron, coherent and three tagged specs, where the
`li6-convolution` BASE of the cell is itself refused, so there is no run to
vary; the matrix says so by name rather than passing silently.

The tables above are `phase_A_numbers.md`'s own; the pinned subset lives in
`python/tests/test_b1_model.py` P10 (`R_WEIGHT_STATUS_QUO`,
`R_WEIGHT_SHARED_R1998`, `R_WEIGHT_NUMERATOR_ONLY_PCT`,
`R_COS2PHI_SHARED_PCT`) and in `tests/test_b1_nuclear.cpp` T18, so a silent
drift fails a suite rather than a reading.

**Option (ii) has no flag and is not getting one.** It is built by hand in
both test files (`_numerator_only_r1998_kernel`), because (ii) is the option
the registry did not take and giving it a flag would be taking a decision.

### A1.9 Counts this section moved

Adding the axis moved these published counts, all by construction (measured after the change, on this machine with `env.sh` sourced):

| count | before (baseline above) | after A1 | why |
|---|---|---|---|
| `knob_provenance` rows on a default run | 66 | **67** | the `r_source` row |
| `knob_provenance` matrix variants | 71 | **74** | three `r_source` cells |
| matrix cells | 547 | **568** | 3 variants × the 7 full specs |
| `test_knob_provenance.py` collected | 582 | **603** | the same 21 |
| refused axes on the CLI's default run | 7 | **8** | `r_source` is refused on `miller` |
| `build/lipolgen_tests` cases / assertions | 401 / 17 240 286 | **402 / 17 240 344** | T18 and its 58 assertions |
| `pytest python/tests` | 926 passed / 112 skipped | **941 / 124** | P10's 6, plus the matrix's 9 + 12 |
| docs gate (strict) | 1206 refs / 96 ranges / 7 external / 0 broken / 6 allowed | **1216 / 96 / 7 / 0 / 6** | the `RSource` row's 10 citations |
| SPDX | 101/101 | **101/101** | no new source file |

No pre-existing C++ assertion, pytest or citation moved: every delta above is
an addition and its count is stated. The nine `validation/reference/*.json`
are byte-identical (`git diff --stat validation/reference/` is empty).

The live statements of them (`include/lipolgen/pipeline.hpp`,
`docs/USAGE.md` §7c and its verbatim banner block, `docs/CONVENTIONS.md`,
`README.md`, `docs/PHYSICS_CHANNELS.md`) are updated. The dated run records
that quote 547/582 (`../run_2026-09-03/SUMMARY.md`,
`../run_2026-09-03/phase_D_numbers.md`) are **left as measured on their own
date** — the tree's rule for a dated record (`phase_E_numbers.md` §E1.5).

---

## A2 — registry row 17: the RC band's LOW-x ANCHOR, `RC_DELTA_LOW_X` = 0.30

**What this section prices, and what it does not decide.** `STATUS.md`
decision row 17 (`../run_2026-09-03/AUTHOR_DECISIONS.md` §B13) asks whether the
tensor-RC band's low-x anchor stays at δ = 0.30, hung at `RC_X_LOW` = 0.01.
Phase B established what the 0.30 rests on and the registry recorded the two
alternatives **with no band re-evaluation** — "unpriced" in its own index. They
are re-evaluated here, by **running the band**, not by scaling one row.
**Nothing is decided: no default moved, no reference JSON moved, and the knob
that does the work (`--rc-delta-low-x`) already existed — this section adds no
option.**

### A2.1 What the 0.30 rests on, and which way it errs

Gakh–Shekhovtsova (hep-ph/0403262, **zero INSPIRE citations**) say the
correction "changes from 10 % to 30 % as compared with the Born contribution"
and scope that to **x ∼ 10⁻³–10⁻²**. Exactly one panel of their Fig. 2 lies in
that window — panel (a), Q² = 0.1 GeV², x = 0.00226–0.00966 — and its two ends
are the 10 % and the 30 %: |δ| = **0.2662** at x = 0.00226 and **0.1133** at
x = 0.00966 (`../run_2026-09-03/phase_B_numbers.md` §B6.3, measured off the
arXiv figure arrays). `RC_X_LOW` = 0.01 sits **above** the panel's top. So:

> **0.113 is the value the panel actually reads at the x nearest the anchor.**
> The shipped 0.30 is **×2.65** that reading and **×1.13** the 0.266 at the
> panel's own bottom. It errs **wide** — the safe direction for a half-width
> and the wrong one for a quoted precision.

The two alternatives are those two readings. This section runs the registry's
own rounded values, **0.266 and 0.113**; the raw arrays read 0.2662 and 0.1133,
so a half-width built on the raw readings would be **+0.075 %** and **+0.27 %**
above the rows below — far under any digit that matters here, and stated so
that nobody re-derives a discrepancy from it.

### A2.2 The three-row table

The configuration is §9's own, so the rows are comparable with everything else
in that section: ⁶Li, `--config 1` (s per nucleon = 3980 GeV²),
`--channel inclusive`, `--plan tensor-thirds --pzz 0.6`, θ_S = 0, defaults
(`ho`, `t-peak`, `a_transfer_frac` 0), **Q² = 5 GeV², P_zz = +1**. The
half-width on `A_zz` is `RcModel::delta(x)·|A_zz|` — the model's own single
call site of `rc_delta`, so the `a_transfer_frac` term applies exactly once.

`A_zz` is **this generator's own ⁶Li b₁** (`Li6B1(MillerB1)`), i.e. the value
corrected in `../run_2026-09-03/phase_B_numbers.md` §B3.2, **not** the
×3.253983 deuteron-b₁ column that `OPEN_ITEMS_SOLUTIONS.md` §9 and
`../run_2026-09-02/phase_C_numbers.md` §8.2 still print (both now flagged in
place). Re-measured here at all four x:

| x | y = Q²/(x s) | `A_zz`(Born), ⁶Li b₁ |
|---|---|---|
| 0.010 | 1.256281e−01 | **−4.528241648849e−04** |
| 0.063 | 1.994097e−02 | **−1.316774298646e−03** |
| 0.100 | 1.256281e−02 | **−1.528923484258e−03** |
| 0.160 | 7.851759e−03 | **−1.280460414420e−03** |

The x = 0.010 and x = 0.100 entries reproduce §B3.2 to every printed digit;
0.063 and 0.160 are new here. δ(x) and the half-width:

| δ_low | δ(0.010) | δ(0.063) | δ(0.100) | δ(0.160) |
|---|---|---|---|---|
| **0.30 — shipped** | 0.30000000 | 0.11080618 | 0.06331262 | 0.01500000 |
| 0.266 — panel bottom (x = 0.00226) | 0.26600000 | 0.09937667 | 0.05754901 | 0.01500000 |
| **0.113 — panel AT the nearest x (0.00966)** | 0.11300000 | 0.04794388 | 0.03161276 | 0.01500000 |

| **band half-width on `A_zz`** | x = 0.010 | x = 0.063 | x = 0.100 | x = 0.160 |
|---|---|---|---|---|
| **(i) 0.30 — shipped** | **1.358472e−04** | **1.459067e−04** | **9.680016e−05** | **1.920691e−05** |
| (iii) 0.266 | 1.204512e−04 | 1.308566e−04 | 8.798804e−05 | 1.920691e−05 |
| **(ii) 0.113** | **5.116913e−05** | **6.313127e−05** | **4.833349e−05** | **1.920691e−05** |
| ratio (iii) ÷ (i) | ×0.886667 | ×0.896851 | ×0.908966 | **×1** |
| ratio (ii) ÷ (i) | ×0.376667 | ×0.432682 | ×0.499312 | **×1** |

(The registry's option letters: **(ii)** is 0.113 and **(iii)** is 0.266,
§B13(b). The rows are ordered by magnitude here, so read the letters and not
the position.)

Three readings, all measured:

1. **The peak does not move — it sharpens.** The band peaks at **x = 0.063 on
   every anchor**, not at the lowest x, because δ(x) is still rising there
   while |A_zz| has not yet fallen. The 0.063 : 0.010 ratio rises
   **1.074 → 1.086 → 1.234** as the anchor falls. So the "largest RC
   systematic sits in the middle of the low-x range" reading survives the
   anchor decision — and gets *more* true under (ii).
2. **The ordering of the RC budget survives every option.** At x = 0.01,
   Q² = 5 the whole radiative tail is ΔA_zz = **−1.769797e−07** (re-measured
   here from `tail_ratio_at`, reproducing §B3.2 exactly), so the band leads it
   by **×767.6 / ×680.6 / ×289.1**. At the other three x the lead is
   ×17 020–×45 500. Even the smallest anchor leaves the band the dominant RC
   systematic by two and a half orders of magnitude.
3. **The high anchor pins x ≥ 0.16 exactly.** 1.920691e−05 at all three, to
   every bit. The entire low-anchor decision lives **below x = 0.16**.

### A2.3 The registry's own arithmetic was right at one x and wrong as a rule

§B13(d) estimated "≈ ×0.89 under (iii) and ≈ ×0.38 under (ii)", calling it
"arithmetic on the band's linear form, not a run". Measured:

* **At x = 0.01 it is exact** — ×0.886667 and ×0.376667, the ratio of the
  anchors themselves, because `rc_delta` returns δ_low **exactly** at
  x ≤ `RC_X_LOW` (checked at x = 0.005 too).
* **Between the anchors it overstates the reduction.** δ(x) is log-linear in
  x between two *fixed* anchors, so moving only the low one moves δ(x) by a
  fraction of the anchor's own move: **×0.4327 at x = 0.063 and ×0.4993 at
  x = 0.100 under (ii)**, against the ×0.3767 at the anchor. The ratio is
  strictly between δ_low/0.30 and 1 at every x in (0.01, 0.16) — gated.
* **At x ≥ 0.16 it is zero, not "≈".**

So "the band scales with δ_low" is a statement about x = 0.01 and nowhere
else, and the estimate that stood in the registry was a lower bound on the
surviving half-width across most of the band's range.

### A2.4 What actually moves in a file

⁶Li inclusive, `--config 1`, 2000 events, seed 99, `--rc tensor-band`,
compared column by column as `ndarray.tobytes()` plus a key-by-key `meta` diff
(the §B6.5 method):

| | (iii) 0.266 vs shipped | (ii) 0.113 vs shipped |
|---|---|---|
| columns that differ, of **51** | `rc_tensor_lo`, `rc_tensor_hi` — **and no others** | the same two |
| `meta` keys that differ, of **59** | `rc_delta_low_x`, `knob_provenance` (which **records** it) | the same two |
| mean \|w_hi − 1\| | 6.153874e−05 → **5.484602e−05** (×0.891244) | → **2.472882e−05** (×0.401842) |
| per-event (w_hi − 1) ratio, min / max / mean | 0.8866667 / 1.0000000 / 0.8931602 | 0.3766667 / 1.0000000 / 0.4123810 |

The per-event ratio's **minimum is the anchor ratio** (events at x ≤ 0.01) and
its **maximum is exactly 1** (events at x ≥ 0.16): the same statement as
§A2.3, seen event by event. No kinematic column, no `weight`, no `rc_tail`
moves — the anchor is a band knob and nothing else.

### A2.5 The tagged band, where the clamp bites: the clipped fraction does not move

The one place the band is not a smooth rescaling is the tagged channels, where
`tau_tag = 1 − n̄/n_M` diverges at the nodes of the M-dependent spectator
density and `RcOptions::band_tau_max` = 1 clamps it. 20 000 events, `--config 1`,
`--pzz 0.6`, ⁶Li at **seed 1** and ⁷Li at **seed 11**.

> **RE-MEASURED 2026-09-15 — the ⁶Li column of this table was measured at
> `cdd8591`, BEFORE the tagged S–D interference fix `a7b3d18`, and does not
> reproduce on the committed tree.** The fix moves the spin-1 (M, cos θ_k)
> draw that `RcModel::tagged_tau` reads, so the ⁶Li tagged clipped count moves
> with it; ⁷Li is J = 3/2 and does **not** move. The pre-fix ⁶Li figures were
> **124 = 0.6200 %** and mean \|w_hi − 1\| 2.587195e−02 / 2.297940e−02 /
> 9.962921e−03; they are superseded by the row below and must not be quoted.
> Every ⁷Li figure reproduces to every printed digit, as do all four
> band-edge rows and the structural claim this section exists for.

| | δ_low = 0.30 | 0.266 | 0.113 |
|---|---|---|---|
| ⁶Li clipped events / 20 000 (re-measured 2026-09-15) | **245 = 1.2250 %** | **245 = 1.2250 %** | **245 = 1.2250 %** |
| ⁷Li clipped events / 20 000 | **520 = 2.6000 %** | **520 = 2.6000 %** | **520 = 2.6000 %** |
| non-band columns vs the 0.30 run | — | **all byte-identical** | **all byte-identical** |
| min `rc_tensor_hi`, both isotopes, `== 1 − δ_low` exactly | **0.700** | **0.734** | **0.887** |
| max `rc_tensor_lo`, both isotopes, `== 1 + δ_low` exactly | **1.300** | **1.266** | **1.113** |
| mean \|w_hi − 1\|, ⁶Li (re-measured 2026-09-15) | 3.384715e−02 | 3.006184e−02 | 1.302792e−02 |
| mean \|w_hi − 1\|, ⁷Li | 8.127222e−02 | 7.218497e−02 | 3.129231e−02 |

**Identical, not merely close** — and the reason is structural, not
statistical: `RcModel::clamp_tau` (`src/core/rc.cpp`) clips |τ| against
`band_tau_max` and **δ never enters it**; δ multiplies the already-clamped τ
afterwards. So the anchor cannot move *which* events clip. What it moves is
the **width of the band on the events that do**: |τ| is clamped to exactly 1
there, so their edges are exactly 1 ∓ δ_low, and the extrapolated band on the
worst-behaved 1.23 % of a tagged run is **[0.700, 1.300]** on the shipped
anchor against **[0.887, 1.113]** on the panel's own reading. That is the real
cost of the anchor at the clip: a factor 2.65 on the width of the least
defensible part of the band.

**Two bookkeeping facts from the same runs, recorded because they are quoted
elsewhere.** The ⁶Li **1.2250 %** is **seed 1** and sits mid-scatter in a
ten-seed scan — 1.0950 / 1.1150 / 1.2000 / 1.2250 / 1.2450 / 1.2800 / 1.2800 /
1.3350 / 1.3450 / 1.3650 % at seeds 4 / 2 / 20260713 / 1 / 1234 / 4242 / 99 /
11 / 3 / 42, **mean 1.2485 %, sd 0.0923** (re-measured 2026-09-15; the
pre-fix scan was 0.4850–0.6200 %, mean 0.5280 %, sd 0.0411, with seed 1 at its
**top** — that reading of seed 1 does not survive the fix either). So it
should be read as one draw from a scatter, not as *the* fraction. The ⁷Li
**2.6 %** is **seed 11**; seed 1 there gives **2.51 %** — both reproduce
unchanged on the fixed build.

### A2.6 The HIGH anchor: the record names no alternative

The brief asks whether the E12-13-011 anchor has one. **It does not, and this
run does not invent one.** `RC_DELTA_HIGH_X` = 0.015 at `RC_X_HIGH` = 0.16 is
the proposal's *"the unpolarized corrections are known to better than 1.5 %"*.
Grepped across `include/`, `src/`, `python/`, `tests/` and `docs/`: there is
**no second reading of it, no `_OPTIMISTIC` partner** (the low anchor has one,
`RC_DELTA_LOW_X_OPTIMISTIC` = 0.19) **and no band edge** anywhere. What the
record does carry is a *provenance* flag and an *action*, neither of which is a
number:

* the 1.5 % lives in the **unpublished proposal only**; the published
  companion (Poudel *et al.*, arXiv:2506.04506, EPJ A **61** (2025)) contains
  no radiative-correction discussion at all — the string "radiative" does not
  occur in it (`../run_2026-09-02/design_C_tensor_rc.md`);
* the only alternative any document states is **"cite it by page for the
  1.5 %, or drop the anchor"** (`design_C_tensor_rc.md`; `CONVENTIONS.md` (a)).

So the x ≥ 0.16 column of §A2.2 is a single number by necessity, and the honest
statement is that the high anchor is **unbanded because there is nothing to
band it against** — not because it is better known than the low one.

### A2.7 What it would cost to adopt (ii) or (iii)

Nothing numerical that is shipped: the band is opt-in (`--rc tensor-band`),
`--rc off` is byte-identical, and **no `validation/reference/*.json` carries
an `rc_*` field at all**, so none of the nine rtol-1e−12 references can move
under any choice. No existing assertion breaks either — every gate on the
anchor is symbolic (`rc_delta(RC_X_LOW) == RC_DELTA_LOW_X`,
`|w − 1| ≤ RC_DELTA_LOW_X`, `worst_lo ≈ 1 − RC_DELTA_LOW_X`), and the two
literal `0.30`s in the suites are a bad-anchor `CHECK_THROWS` argument and an
explicitly-passed `rc_delta_low_x=0.30` construction, neither of which pins the
default. The cost is **the constant, the pinned table, and the sites that
state the claim**: `include/lipolgen/rc.hpp` (`RC_DELTA_LOW_X`, `RC_X_LOW` and
the honest-flags block), `python/tests/test_rc_low_x_anchor.py`'s
`HALF_WIDTH`, the `--rc-delta-low-x` help and the banner reason,
`docs/USAGE.md` §7b, `docs/OPEN_ITEMS_SOLUTIONS.md` §9,
`docs/PHYSICS_CHANNELS.md`, `docs/CONVENTIONS.md` (b) and
`design_C_tensor_rc.md` Q7.

### A2.8 How to reproduce

```
cd LiPolGen && source env.sh && cmake --build build -j
python -m pytest python/tests/test_rc_low_x_anchor.py -q      # 5 passed
```

The §A2.2 tables, from the shipped Python API alone (the §B6.2 recipe with the
anchor varied instead of `a_transfer_frac`):

```python
import lipolgen as lg
def mk(**kw):
    cfg = lg.make_config(isotope="6Li", channel="inclusive", config=1,
                         events=2000, seed=99, rc="tensor-band", **kw)
    return lg.Pipeline(cfg, lg.tensor_thirds_plan(0.0, 0.6))
p0 = mk(); k, s, m0 = p0.dis_sampler.kernel, p0.dis_sampler.s, p0.rc_model
for x in (0.01, 0.063, 0.10, 0.16):
    t = k.tables(x, 5.0)
    azz = lg.azz(t.b1, t.f1, t.f2, x, 5.0 / (x * s), t.b2, 0.0)   # 6Li b1
    rU = m0.tail_ratio_at(x, 5.0, 0.0)
    rT = m0.tail_ratio_at(x, 5.0, 1.0) - rU                       # the tail
    print(x, azz, (azz + 2 * rT) / (1 + rU) - azz,
          [(mk(rc_delta_low_x=d).rc_model.delta(x) * abs(azz))
           for d in (0.30, 0.266, 0.113)])
```

and §A2.5, with `channel="tagged-alpha"`, `events=20000`, `seed=1` (⁶Li) —
`meta["rc_clipped_band_events"]` and `cols["rc_tensor_hi"].min()`. The ⁷Li row
needs the API plan of `test_rc.py`'s J = 3/2 tagged case, at `seed=11`.

### A2.9 The default did not move — measured, not argued

This section changed no code path: what it touched in `include/lipolgen/rc.hpp`
is comment only, in `python/lipolgen/cli.py` two `help=` strings and one banner
line, and in `src/core/pipeline.cpp` two `band_row` **reason strings**. That
was checked the §B6.5 way rather than asserted — the two `band_row` strings
were reverted in place, the library rebuilt, and a 2000-event ⁶Li inclusive run
at seed 4242 regenerated on **both** `--rc off` and `--rc tensor-band` and
compared as `ndarray.tobytes()` per column plus value by value on `meta`:

| compared | entries | differing |
|---|---|---|
| generated columns (48 at `--rc off` + 51 at `--rc tensor-band`) | **99** | **0** |
| `meta` values (32 at `--rc off` + 59 at `--rc tensor-band`) | **91** | **1** |

The one difference is `meta["knob_provenance"]` on the `--rc tensor-band` run
— the **reason text** of the `rc_delta_low_x` and `rc_delta_high_x` rows, i.e.
the sentences this section wrote. Every numeric `meta` value is identical, and
`--rc off` is identical in **both** columns and `meta` (at `rc = off` those two
rows carry the "no `RcModel` is built" reason, which was not touched). The
files were then restored byte-for-byte and rebuilt.

### A2.10 Counts this section moved

| count | before (after A1 + A3) | after A2 | why |
|---|---|---|---|
| `build/lipolgen_tests` cases / assertions | 404 / 17 240 384 | **404 / 17 240 384** | unchanged — no C++ test changed, only a banner string |
| `pytest python/tests` | 987 passed / 127 skipped | **992 / 127** | `test_rc_low_x_anchor.py`, 5 tests |
| docs gate (strict) | 1220 refs / 95 ranges / 7 external / 0 broken / 6 allowed | **1220 / 95 / 7 / 0 / 6** | unchanged — **no citation added**; the comment edits in `rc.hpp`, `cli.py` and `pipeline.cpp` moved line numbers, so `--fix` re-anchored every reference it could resolve, three **use-site** citations that `--fix` cannot relocate (their symbols are declared in no file) were bumped by hand and verified against their recorded line text — `cli.py:1326 → 1343` (`rc_model`), `cli.py:1092 → 1109` (`set_pythia_hadronizer`), `pipeline.cpp:3548 → 3557` (`hadronizer`) — and `--record-ranges` re-fingerprinted the two blocks this section edited in place, `rc.hpp:282-361 → 326-408` and `USAGE.md:2029-2182 → 2029-2204`. The recorder reported **0 changed** on every pre-existing fingerprint |
| SPDX | 101/101 | **101/101** | `python/tests/*.py` is outside the covered globs (the new file carries the header anyway) |
| `knob_provenance` rows on a default run | 68 | **68** | **no new option** — `--rc-delta-low-x` was already registered |

No default moved, no reference JSON moved
(`git diff --stat validation/reference/` empty), and no pre-existing assertion
or citation was edited.

---

## A3 — registry row 20 / §15.5 **D13**: honouring `--pzz` on `helicity-flip`

**What this section prices, and what it does not decide.** `STATUS.md` decision
row 20 (`../run_2026-09-03/AUTHOR_DECISIONS.md` §B18) asks what to do about
`--pzz` on the `helicity-flip` plan, and offers three answers: **(i)** leave it
documented (applied 2026-09-05: the banner prints the fill's own moments);
**(ii)** refuse `--pzz` with that plan; **(iii)** honour it. The run recorded
that (iii) "would move a shipped fill" and that **its effect on any observable
was not measured**. (iii) is now implemented, **opt-in, default unchanged**,
and that effect is measured here. **Nothing is decided: no default moved, and
all nine `validation/reference/*.json` are byte-identical.**

### A3.1 What was wrong, exactly

`--plan helicity-flip --pzz 0.6` at `--pz 0.7` builds
`helicity_flip_plan(j, pz, pe)` with `HelicityFlipOptions` untouched, so
`use_explicit_pzz` stays false and the fill is `populations_maxent(j, pz)` —
the max-entropy ladder. The typed 0.6 reaches nothing. At J = 3/2 the ladder
fills **T = 0.4**; at J = 1 it fills **P_zz = 0.4094026279413209**. Since
2026-09-05 the banner has printed the fill's own moments, so the substitution
is not silent, but no switch honoured the typed value.

0.6 at `--pz 0.7` is outside the J = 3/2 domain **by 0.02**: the four
populations are fixed uniquely by (1, P_z, T, R₃), so positivity is
1.8|P_z| − 1 ≤ T ≤ 1 − 0.6|P_z| at R₃ = 0, i.e. **0.26 ≤ T ≤ 0.58** here, and
at T = 0.58 exactly p(−½) = 0. It is **inside** the spin-1 domain
3|P_z| − 2 ≤ P_zz ≤ 1, i.e. **0.1 ≤ P_zz ≤ 1**, so the same command is
physical on ⁶Li and not on ⁷Li.

### A3.2 The switch, and why the default cannot move

`--pzz-mode {ladder,typed}`, default `ladder`; Python
`make_plan(..., pzz_mode="ladder")`; `PZZ_MODES` in `lipolgen/__init__.py`.
It is **the name of `HelicityFlipOptions::use_explicit_pzz`'s two branches and
nothing else** — `make_plan` is the one place the string becomes the bool, and
`ladder` does not touch the options object at all, so the default takes the
branch that was already there. That is bit for bit **by construction**, not by
an argument that two code paths agree, and `test_pzz_mode.py` P1 checks it as
identity of doubles (populations, `pzz_true`, `pz_true`) and as byte equality
of every generated column of a 400-event run against the pre-flag call
`_l.helicity_flip_plan(j, pz, pe)`.

`typed` **refuses** a value outside the plan's domain rather than clamping to
its edge — a clamped fill would publish an alignment nobody typed. The refusal
is `spin32_populations`' (which has named its edges since 2026-09-05) and,
newly, `spin1_populations`', which until now threw the useless "unphysical
(pz, pzz): negative population" and now names the offending m, its negative
population and both edges. The CLI turns either into a `SystemExit` carrying
that message instead of a traceback.

### A3.3 The standard configuration

Inclusive channel, beam config 1, **100 000 events, seed 20260713**,
`--plan helicity-flip --pz 0.7 --pe 0.7`, the shipped `--b1-model miller`,
`--pzz 0.5` — a value **inside both domains**, so the two fills can be
compared at all. `sigma_pb` is the cell-grid integral and does not depend on
the event count (checked at 200 / 2 000 / 100 000: the same doubles), which is
why `test_pzz_mode.py` P4 pins it at 200 events.

| | ⁶Li, J = 1 | ⁷Li, J = 3/2 |
|---|---|---|
| ladder populations (m = +J … −J), Python `repr` | 0.7515671046568897, 0.1968657906862264, 0.05156710465688399 | 0.6750000000000008, 0.2249999999999996, 0.07499999999999965, 0.02499999999999981 |
| typed populations | 0.7666666666666666, 0.16666666666666666, 0.06666666666666672 | 0.69, 0.2300000000000001, 0.019999999999999934, 0.060000000000000074 |
| recorded alignment, ladder → typed | 0.4094026279413209 → 0.5 (**+22.129162 %**) | 0.40000000000000147 → 0.49999999999999994 (**+25 %**) |

### A3.4 Every observable that reads the fill

Ladder → typed, at the configuration of A3.3.

| observable | ⁶Li (J = 1) | ⁷Li (J = 3/2) |
|---|---|---|
| σ_pb | 591783.2520093301 → 591769.3290332475 pb, **−0.0023527 %** | 590952.426415097 → 590952.4264150971 pb, **+1.970 × 10⁻¹⁶ relative (1 ulp)**; not bit-identical |
| σ[apar+] | 591443.485190121 → 591429.562214039, **−0.0023541 %** | 591105.66165053 → 591105.66165053, **3.939 × 10⁻¹⁶ relative (2 ulp)** |
| σ[apar−] | 592123.018828539 → 592109.095852456, **−0.0023514 %** | 590799.191179664 → 590799.191179664, **exactly equal** |
| counts per category | 49972 / 50028 → unchanged | 50013 / 49987 → unchanged |
| A_∥ from the rate, (σ₊−σ₋)/(σ₊+σ₋)/(P_e P_z) | −0.00117171560617921 → −0.00117174317396206, i.e. **−2.757 × 10⁻⁸ = −4.2717 × 10⁻⁶ σ_stat** | +4.019 × 10⁻¹⁶ absolute = **+6.23 × 10⁻¹⁴ σ_stat** |
| δ(A_∥) = `err_a_parallel(N, P_e, P_z)` | 0.00645362787789465, **does not move** (it has no P_zz) | same, does not move |
| δ(A_zz), δ(cos 2φ) at fixed N (∝ 1/P_zz) | **−18.119474 %** | **−20.000000 %** |
| generated columns that move | **15 of 48** (`cell, e_prime, eta_e, kp, m, m_ion, nu, phi, pzz, q2, struck_pdg, theta_e, w2, x, y`; corrected 2026-09-06 — `R`, `m_struck` and `xL` are all-NaN, byte-identical columns on the inclusive channel that a NaN-blind `!=` counted as moved; the byte comparison `test_pzz_mode.py` P1 and the knob matrix use gives 15) | the same 15 |
| events landing in a different (x, Q²) cell | **1201 of 100 000 = 1.2010 %** | **106 of 100 000 = 0.1060 %** |
| mean Q² of the accepted sample | 4.3162908488542 → 4.29061031437528 GeV², **−0.594968 %** | 4.21491916332097 → 4.2235038104617, **+0.203673 %** — *all of it from the 106 re-celled events, see A3.5* |

**A 1σ shift of A_∥ would need N ≈ 5.5 × 10¹⁵ events on ⁶Li** and
**≈ 2.6 × 10³¹ on ⁷Li** (δ ∝ 1/√N applied to the shifts above). Neither is an
experiment.

### A3.5 Why ⁷Li moves nothing, and why it still moves 106 events

⁷Li's rank-2 sector is identically zero (`rank2_input_report`; there is no
`b1_32`), so the only fill moments `InclusiveKernel::amplitudes` sums are
Σ p_m and ⟨J_z⟩/J — and **both fills honour `--pz`**, which is the whole
content of the result. Measured from the populations themselves:

| fill | Σ p_m | ⟨J_z⟩/J | T |
|---|---|---|---|
| ladder | 0.99999999999999978 | 0.70000000000000107 | 0.40000000000000147 |
| typed 0.5 | 1 | 0.69999999999999984 | 0.49999999999999994 |

The two vector moments differ by **1.2 × 10⁻¹⁵**, and that is the entire
difference the kernel can see. T differs by 0.1 and multiplies zero.

**It is nevertheless not a no-op on the file.** Cell selection is a discrete
function of weights that differ in the last bit, so 106 of 100 000 events land
in another (x, Q²) cell and 15 of the 47 columns change (byte comparison; three all-NaN columns do not). A ⁷Li reference built
under one mode would not reproduce byte for byte under the other, even though
no physical number moved. That is the honest form of "(iii) moves a shipped
fill" on ⁷Li: it moves the **sample**, not the physics.

**And the re-celling is the WHOLE of the sample-level difference, measured.**
Every event whose `cell` label is unchanged carries a **bit-identical** `x`
and `q2` under both modes (`np.array_equal` over the 99 894 of them; the same
holds for the 98 799 on ⁶Li). The mean-Q² shift of A3.4 is therefore entirely
the 106 movers: Σ Q² over them goes 3546.44 → 4404.91, a change of **858.465
GeV²**, and the change of the total Σ Q² over all 100 000 events is
**858.465** — the same number to every digit. (On ⁶Li: −2568.05 over the 1201
movers and −2568.05 over the whole sample.) A single event crossing into the
high-Q² tail moves the mean visibly because the window reaches Q² = 1030 GeV²;
that is what a 0.2 % mean shift with 106 movers means, and it is not a
physics difference.

### A3.6 On ⁶Li the size is the b₁ model's, not the window's

`li6-convolution` needs `--x-max 0.95`, so the window has to be controlled
before the model can be blamed. Same runs, same `--pzz 0.5`:

| configuration | σ ladder → typed | shift | A_∥ shift | cells moved |
|---|---|---|---|---|
| `miller`, `x_max` 1.0 (shipped) | 591783.252009 → 591769.329033 | **−0.0023527 %** | −4.2717 × 10⁻⁶ σ | 1201 |
| `miller`, `x_max` 0.95 | 591783.244679 → 591769.321734 | **−0.0023527 %** | −4.2717 × 10⁻⁶ σ | 1201 |
| `li6-convolution`, `x_max` 0.95 | 591843.294739 → 591842.660370 | **−0.00010719 %** | −1.9459 × 10⁻⁷ σ | 1290 |
| `cdks`, `x_max` 0.95 | 591844.391798 → 591844.000199 | **−0.00006617 %** | −1.2012 × 10⁻⁷ σ | 1291 |

The window changes the shift in the **8th** significant figure; the b₁ model
changes it by a factor **22** (`li6-convolution`) to **36** (`cdks`) below
`miller`. So "what honouring `--pzz` costs the ⁶Li rate" is a statement about
the b₁ backend, and `miller` — the shipped default — is the configuration
where it costs most.

### A3.7 The refusal, shown

The registry's own line, `--plan helicity-flip --pzz 0.6` at `--pz 0.7`:

```
$ lipolgen-run --isotope 7Li --plan helicity-flip --pzz-mode typed --events 10
spin32_populations: unphysical (pz = 0.7, t = 0.6, o = 0) -- p(m = -1/2) =
-0.005 < 0.  The four populations are fixed UNIQUELY by (1, pz, t, o), so this
is a DOMAIN, not a solver failure: ...  At o = 0 that is 0.26 <= t <= 0.58 at
this pz (so |pz| <= 5/6 = 0.833333 at all).  The J = 3/2 domain is SMALLER
than the spin-1 one -- (pz, pzz) = (0.7, 0.6) is inside it for spin 1 and
outside here by 0.02 in t
```

and the same 0.6 **runs on ⁶Li**, where it is inside the domain: σ
591783.2520093301 → 591753.9610642146 pb = **−0.0049496 %**, A_∥
**−8.9869 × 10⁻⁶ σ_stat**, **2541 of 100 000 events (2.541 %)** re-celled,
alignment 0.409403 → 0.6. At the default `--pzz-mode ladder` both isotopes
still fill the ladder and neither refuses.

The spin-1 refusal now names its own edges too, which it did not before this
section (checked against the function on a 41 × 41 (P_z, P_zz) grid, **0
disagreements**):

```
spin1_populations: unphysical (pz = 0.7, pzz = -0.5) -- p(m = -1) = -0.1 < 0.
The three populations are fixed UNIQUELY by (1, pz, pzz), so this is a DOMAIN,
not a solver failure: p(0) = (1 - pzz)/3 and p(+-1) = (2 +- 3 pz + pzz)/6,
whence 0.1 <= pzz <= 1 at this pz (so |pz| <= 1 at all).  The J = 3/2 domain
is SMALLER: at J = 3/2 the same pz admits 0.26 <= T <= 0.58
```

### A3.8 The provenance rows, and what they are measured from

`--pzz-mode` gets a `knob_provenance` row whose status is a **measurement**,
not a table of plan names: the row rebuilds the *other* mode's fill at this
run's own (J, P_z, P_e) and counts the populations that differ as doubles —
the `RouteReach` rule applied to a fill. `read` under `helicity-flip`,
`not-read` (labelled `not read by plan tensor-thirds`, …) under the three
tensor factories, which have no ladder branch. This matters at J = 1, where
the populations are fixed uniquely by (P_z, P_zz): a typed value that
reproduced the ladder's fill exactly would have to be reported `not-read`, and
a rule keyed on the plan name would get that wrong.

> **Qualified 2026-09-15 — "a measurement" holds where the other fill can be
> BUILT, and on three cells it cannot.** On the ⁶Li helicity-flip specs the
> reason text does carry the count, verbatim: *"MEASURED here by rebuilding the
> other mode's fill: 3 of 3 populations are different doubles"*. On the **three
> J = 3/2 helicity-flip specs** — `inclusive-7Li/helicity-flip`,
> `tagged-7Li-alpha/helicity-flip` and its `pe = 0` sibling — the typed fill is
> **outside the plan's domain** at the matrix's own (P_z = 0.7, P_zz = 0.6):
> `spin32_populations` refuses it (*"unphysical (pz = 0.7, t = 0.6, o = 0) —
> p(m = −1/2) = −0.005 < 0"*). The row there is `read` **because the axis is
> consulted and the refusal is quoted in the reason**, not because two
> population vectors were compared — no populations were counted on those three
> cells. The reason text has always said so; it is this section's
> generalisation that was stronger than the J = 3/2 cell, and the matrix test
> accepts the arrangement because a refused variant leaves the baseline row
> unjudged (`test_knob_provenance.py:586-597`).

The `--pzz` row moves with it. It was `not-read` under `helicity-flip`
unconditionally; it is now `not read by plan helicity-flip at --pzz-mode
ladder` (with the measurement in its reason: *this run's populations ARE
`populations_maxent(J, --pz)`, double for double*) and **`read`** at
`--pzz-mode typed`. Whether the fill came from the ladder is likewise
measured, from the plan, so a C++ caller that supplies no `pzz_mode` is
classified rather than guessed at.

Both rows are in `meta["knob_provenance"]` and in the banner's provenance
block, and the banner's `fill` line now names the mode and the fill the other
mode would have built:

```
  fill helicity-flip: J = 1.5, P_z = 0.7, T = 0.4
    --pzz-mode ladder (the default): the max-entropy fill at --pz, so
    --pzz = 0.6 is NOT read -- --pzz-mode typed would build the fill at
    T = 0.6 instead (refused, with the edge named, where that is outside the
    plan's domain)
```

### A3.9 What (iii) costs, stated for the decision

* **On ⁷Li — the isotope row 20 is about — (iii) buys the recorded alignment
  and costs no observable.** σ to 1 ulp, A_∥ to 6 × 10⁻¹⁴ of its own error.
  What it does cost is byte-reproducibility: 0.106 % of the sample re-cells.
* **On ⁶Li (iii) does move the rate**, by −0.0023527 % on the shipped b₁
  model and by 22–36× less on the opt-in ones. For scale, measured at the same
  window (`x_max` 0.95, ladder fill): the b₁ MODEL itself moves σ by
  **+0.0101473 %** (miller → li6-convolution, 591783.244679 → 591843.294739)
  and **+0.0103327 %** (miller → cdks, → 591844.391798), so the fill choice is
  **0.23 of the b₁-model spread** on this observable — the same order, not a
  negligible one. Against the statistics it is nothing: 4 × 10⁻⁶ of A_∥'s own
  error at the standard 100 000 events.
* **What (iii) really changes is the DIVISOR.** The alignment the file records
  is what every tensor estimator divides by, and it moves +25 % (⁷Li) /
  +22.13 % (⁶Li) — worth −20 % / −18.12 % on δ(A_zz) and δ(cos 2φ) at fixed
  N. A programme that wants a stated P_zz on an A_∥ fill gets it here; nothing
  else does.
* **And (iii) makes a working command fail.** `--plan helicity-flip --pzz 0.6`
  at `--pz 0.7` runs today on ⁷Li (filling 0.4) and is REFUSED under
  `--pzz-mode typed`. That is the intended behaviour — the alternative is
  clamping, which publishes an alignment nobody typed — but it is the one
  user-visible cost of ever making `typed` the default.
* **Which references would move if it became the default**: one of the nine,
  `validation/reference/bookkeeping.json`, which pins five `helicity_flip_plan`
  fills at rtol 1e-12 — four from the ladder branch
  (`helicity_flip_j12`, `helicity_flip_j1_maxent_anchor`,
  `helicity_flip_j32_maxent_anchor`, `helicity_flip_j1_tilted_offset`) and one
  already from `use_explicit_pzz = true` (`helicity_flip_j1_explicit_pzz`).
  The other eight carry no fill. **At the opt-in default all nine are
  byte-identical.**

The row therefore does not close on "no measured effect": on ⁷Li there is none
to five significant figures, and on ⁶Li there is one and it is 2 × 10⁻⁵ of the
rate. What the row now has to weigh is the **alignment** against the
**refusal**.

### A3.10 How to reproduce

```
cd LiPolGen && source env.sh && cmake --build build -j
build/lipolgen_tests -tc="*pzz-mode*"                                 # 1 case, 33 assertions
build/lipolgen_tests -tc="*spin-1 domain*"                            # 1 case,  7 assertions
python -m pytest python/tests/test_pzz_mode.py -q                     # 25 passed
python -m pytest python/tests/test_knob_provenance.py -k "pzz" -q     # 33 passed, 3 skipped
```

The 3 skips are the three J = 3/2 specs on the `--pzz 0.6 -> 0.4 at
--pzz-mode typed` cell, whose BASE is the typed fill at `--pzz 0.6` and is
refused — the domain refusal of A3.7, reached from the matrix and named there
rather than passing silently.

The tables of A3.3–A3.7 are this file's own; the subset that must not drift is
pinned in `python/tests/test_pzz_mode.py` (`PRICE`, the ⁷Li ulp statement, the
1/P_zz error ratios and both refusal messages) and in
`tests/test_bookkeeping.cpp` / `tests/test_spin.cpp`, so a silent drift fails
a suite rather than a reading.

### A3.11 Counts this section moved

| count | before (after A1) | after A3 | why |
|---|---|---|---|
| `knob_provenance` rows on a default run | 67 | **68** | the `pzz_mode` row |
| `knob_provenance` matrix variants | 74 | **76** | `pzz_mode`, and `pzz` under `typed` |
| matrix cells | 568 | **592** | 2 variants × (7 full + 5 plan-axis) specs |
| `test_knob_provenance.py` collected | 603 | **627** | the same 24 |
| refused axes on the CLI's default run | 8 | **8** | `--pzz-mode` is LABELLED, not refused |
| `build/lipolgen_tests` cases / assertions | 402 / 17 240 344 | **404 / 17 240 384** | two cases, 40 assertions |
| `pytest python/tests` | 941 passed / 124 skipped | **987 / 127** | 25 in `test_pzz_mode.py` + 24 matrix cells, 3 of which skip on the refused `typed` base |
| docs gate (strict) | 1216 refs / 96 ranges / 7 external / 0 broken / 6 allowed | **1220 / 95 / 7 / 0 / 6** | four new citations (`--pzz-mode` twice, `PZZ_MODES`, `KnobRunContext::pzz_mode`); the `make_plan` range PHYSICS_CHANNELS row 81 used to cite for "the CLI cannot reach the explicit branch" is dropped, because that sentence is no longer true |
| SPDX | 101/101 | **101/101** | no new file in the covered set |

No pre-existing C++ assertion, pytest or citation moved: every delta above is
an addition and its count is stated. The nine `validation/reference/*.json`
are byte-identical (`git diff --stat validation/reference/` is empty).
