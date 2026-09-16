# Phase D1 — the ⁷Li rank-2 (tensor) input: the zero, the α–t construction, and what can be gated

Task **D1-research** of `docs/open_items/run_2026-09-03/PLAN.md` §Phase D.
**READ-ONLY on the code.** No source file, no data file and no reference JSON was
touched. Every number below was produced by running the shipped library from a
clean tree at commit `903fcc9`; §9 gives the script that reproduces all of them.
Where a number is *not* measured, the sentence says so.

Baseline reproduced after the work, from the same clean tree:

* `build/lipolgen_tests` — **390 test cases, 17217416 assertions, 1 skipped, 0 failed**
* `python -m pytest python/tests -q` — **345 passed**
* `python3 validation/check_physics_channels_links.py` — **1094 references checked (strict), 97 ranges, 6 external, 0 broken, 6 allow-listed**

Companion reading: `docs/theory/SPIN32_FINITE_GAMMA.md` (the J = 3/2
decomposition, the rank-≤2 theorem, the sign of `b1_32`),
`include/lipolgen/b1_nuclear.hpp` (the ⁶Li α–d four-term convolution this note
transposes), `docs/open_items/run_2026-09-02/design_D_b1_li6.md` (its design),
`docs/open_items/vmc_reconciliation.md` (the VMC file inventory).

---

## 0. The headline

**The premise is true. Verified three ways.** A ⁷Li inclusive run's entire
rank-2 sector — the tensor term in the φ-averaged rate, the cos 2φ gluon
transversity amplitude, and therefore A_zz — is **exactly, bit-for-bit zero**,
because `default_inclusive_kernel` fills `b1_32_func / b2_32_func /
delta_32_func` **for no spin**, and the whole `b1_model` machinery is refused
for anything but ⁶Li. It is documented in two internal places and in **no place
a user of the run surface will see**; worse, the npz metadata of such a run
records `b1_model = "miller"`, which is a backend that did not run.

**The honest offline construction exists and is cleaner than ⁶Li's.** ⁷Li's
α + t channel is a **single L = 1 partial wave** (parity plus L ⊗ ½ = 3/2 admits
no other), the triton is J = ½ and the α is J = 0, so **the whole of b₁(⁷Li) in
the two-cluster picture is orbital** — the exact analogue of ⁶Li's terms (2d) +
(2α) with no analogue of its dominant term (1) and no analogue of its term (3).
The alignment coefficient is a pure Clebsch–Gordan number, **exact**, and it is
the *same* coefficient the repository already measures as the ⁷Li polarimeter
⟨P₂(cos θ_k)⟩ = −T/5.

**Measured leading estimate** (Q² = 2.5, shipped `ToyF2` / `r_sigma_lt` / CDKS
Eq. (17) κ = 1, S_αt = 1.0084, the in-tree `momenta/li7_at3.momentum` overlap):

> **b₁(⁷Li)/nucleon = +1.246e−04 / +5.700e−05 / +3.748e−05 / +3.310e−05 /
> +5.485e−05 / +1.040e−04 at x = 0.05 / 0.10 / 0.20 / 0.30 / 0.50 / 0.70**,
> i.e. b₁/F₁ = 3.6e−05 … 1.04e−02, giving |A_T| between 0.004 % and 1.0 %.

> **b₁(⁷Li)_αt / b₁(⁶Li)_orbital = 2.99 ± 0.02** over the whole window and over
> every unpolarised-input, R and κ variation tried — one exception, x = 0.50
> with MSTW, where the ⁶Li denominator nearly vanishes. **⁷Li's tensor input is
> not suppressed. It is three times ⁶Li's orbital term** — because it is not
> paid for by a small D-state amplitude and does not suffer ⁶Li's node
> cancellation. §5.3 derives the 3 as the ratio of the two kernels' **first
> moments about z = 1**, which is why it is F₁-independent.

**A real A = 7 gate exists, and it passes.** There is no A = 3 analogue of the
A = 2 gate — ³H and ³He are J = ½ and have **no rank-2 structure function at
all** — but ⁷Li has a large *measured* rank-2 observable that the same wave
function predicts with no free parameter:

> **Q(⁷Li)_α–t = −(2/5) Z_eff ⟨r²⟩ = −3.4851 fm²** (⟨r²⟩ = 12.5348 fm²,
> r_rms = 3.54 fm, Z_eff = 0.695075), against the measured **−4.00(3) fm² →
> ratio 0.871**, converged to 0.1 % on the grids and ±0.3 % on the VMC MC band.
> The ⁶Li sign gate, from the same code path, gives **−0.3333 fm² against
> −0.0818 fm², a factor 4.08 too large**. **The ⁷Li wave-function input is
> validated by a measured moment an order of magnitude better than ⁶Li's is.**
>
> **COMMITTED 2026-09-06** (run 2026-09-06 §B3, D11 paid). The reference is
> now `LI7_QUADRUPOLE_FM2` = **−4.06 fm²** (`rc.hpp`, TUNL's A = 5, 6, 7
> evaluation — the same one `LI6_QUADRUPOLE_FM2` comes from), so the primary
> ratio is **0.858389** (0.871265 against the −4.00(3) this note quoted), and
> the whole calculation is shipped code (`li7_alpha_t_quadrupole`) that
> reproduces this note's own §9 recipe **bit for bit**, pinned by doctest T13
> and pytest G8. The ⁶Li contrast re-measures as **4.0748**, so "an order of
> magnitude better" is **21.7×**, measured. **It gates the wave function; it
> is not a b₁ verdict and b₁(⁷Li) is still unimplemented.**

**And the recommendation is nevertheless *not* "implement it now".** The
convolution's sign is not predicted by anything in this tree. Swapping the
unpolarised backend `ToyF2 → MSTW2008 LO` — the very swap the A = 2 gate says
you must make before quoting a ⁶Li number — **flips the sign of b₁(⁷Li) at four
of the six x points** and moves it by up to 348 %. A second, independent
in-tree decomposition (⁷Li = ⁶Li + n, `li7_li6_n/li6n_31.table`) gives an
estimate of the opposite sign and up to 8× the magnitude. §8 recommends
**path (c)**: the loud zero and the four run-surface defects **now**, the
quadrupole gate as a committed test **next**, and the b₁ itself only after the
unpolarised-backend decision is made once, for both isotopes.

---

## 1. The finding, verified

### 1.1 Where the zero is

`InclusiveKernel::tables` (`src/core/xsec.cpp:96-117`) dispatches the rank-2
slots on the **kernel's ion spin**: spin 1 reads `b1_func_ / b2_func_ /
delta_func_`, spin 3/2 reads `b1_32_func_ / b2_32_func_ / delta_32_func_`
(`include/lipolgen/xsec.hpp:226` (anchor read 2026-09-15 "b1_32_func, b2_32_func")), and an unset slot yields `0.0`. The
`rank2` branch of `tensor_amplitudes` (`src/core/xsec.cpp:251-278`) is
**entered** for ⁷Li — the geometry `tensor_moments(m)` is non-zero — but every
structure function it multiplies is zero, so the tensor term of `w_avg` and the
whole `a2` are zero.

`default_inclusive_kernel` (`src/core/pipeline.cpp:744-805`) fills the rank-2
slots **only inside `if (std::fabs(ion.spin - 1.0) < 1e-9)`**
(`src/core/pipeline.cpp:690` (as of adec442)). There is no spin-3/2 branch. The Δ (gluon
transversity) slot is set on the same line, so ⁷Li loses cos 2φ as well.

`PipelineConfig::validate` then refuses every `b1_model` other than `Miller` on
⁷Li, with an accurate message (`src/core/pipeline.cpp:592-600`): *"is 6Li ONLY,
got isotope 7Li — spin 3/2 has no rank-2 input here (the 7Li rank-2 slots are
empty by design; there is no published b1 for it)"*. So the flag cannot fill
them either.

### 1.2 Measured at the kernel

```
6Li spin 1.0   f1=0.559097 f2=0.260234 g1=0.0249877 b1=0.00133502 b2=0.000534009 delta=0.00141305
   m=+J, theta_S=0   w_avg=-6.84013e-04  a1=0  a2=0
   transverse        w_avg=+3.42006e-04  a1=0  a2=-1.08596e-03
7Li spin 1.5   f1=0.552621 f2=0.25722  g1=0.0242961 b1=0        b2=0           delta=0
   m=+J, theta_S=0   w_avg=0             a1=0  a2=0
   transverse        w_avg=0             a1=0  a2=0
   tensor_moments(+3/2) = (1.0, 3.0) ;  tensor_moments(1/2) = (-1.0, -3.0)
```

(x = 0.2, Q² = 5, s = 4000 GeV², `default_inclusive_kernel(ion, B1Model::Miller)`.)
The alignment weights Q_NN are correct and non-zero for ⁷Li; only the structure
functions are missing.

### 1.3 Measured at the run surface

A ⁷Li A_zz plan has to be built by hand (§1.4). Two categories, pure |m| = 3/2
(T = +1) against pure |m| = ½ (T = −1), unpolarised beam, θ_S = 0, 60 000
events, seed 11, `--isotope 7Li --channel inclusive`:

```
7Li  sigma per category (pb): [590952.42641509, 590952.42641509]
     A_T = (s+ - s-)/(s+ + s-) = 0.0            <- exactly, both doubles identical
6Li  same construction        : (591692.4891880794, 592153.5282590685)
     asymmetry = -3.8944e-04
```

The two ⁷Li category cross sections are **the same double**. There is no
"small", there is no rounding: the tensor sector is absent, not suppressed.

### 1.4 Is the zero loud anywhere? No.

| Surface | What a ⁷Li user sees |
|---|---|
| `--isotope` help (`python/lipolgen/cli.py:69`) | `"6Li, 7Li, d"` — nothing |
| Run banner (`python/lipolgen/cli.py:561-598`) | Channel, optics, σ per category. The b₁ block prints **only** when `b1_model != Miller`, which ⁷Li can never reach — so **no line at all** |
| npz / HFS `meta` (`python/bindings.cpp:469-471` (as of adec442)) | **`b1_model = "miller"`** — a backend that did not run. `b1_band_scale = 1.0`, `b1_alpha_d_dwave_weight = 1.0`, `b1_unpol = "toy"`, all likewise |
| `PipelineConfig::validate` | Accurate, but only fires if the user *asks* for `--b1-model` (`src/core/pipeline.cpp:592-600`) |
| `docs/PHYSICS_CHANNELS.md:125` | Says it plainly: *"for spin 3/2 no slot is set, so the whole tensor and cos 2φ sector of a ⁷Li inclusive run vanishes."* An internal reference table |
| `docs/USAGE.md:269` (anchor read 2026-09-15 "spin 3/2 and has") | *"`7Li` is spin 3/2 and has no rank-2 input here."* — inside the `--b1-model` bullet, about the **flag**, not about running ⁷Li |

**The answer to the task's question: a user asking for ⁷Li A_zz gets zeros with
no warning, and a metadata block that names a b₁ model.** The `meta` row is the
sharper defect: it is exactly the failure mode `PipelineConfig::validate`
already refuses for `--b1-band-scale` on `miller` — *a knob that did not run
recorded as if it had* — applied here to the model name itself.

### 1.5 Three more run-surface defects found on the way

**F1 — no ⁷Li tensor run plan exists.** `tensor_thirds_plan`,
`transverse_tensor_plan` and `tensor_flip_plan` (`src/core/bookkeeping.cpp:112`,
`:126`, `:134`) all hard-code spin 1. Only `helicity_flip_plan` takes `j`.

**F2 — and `make_plan`'s refusal gives wrong advice.**
`python/lipolgen/__init__.py:477` (as of adec442) refuses `tensor-thirds` at J = 3/2 with
*"needs helicity-flip or transverse-tensor"*. Measured:

```
tensor-thirds      REFUSED: ... J = 1.5 needs helicity-flip or transverse-tensor
transverse-tensor  built OK, categories j = [1.0]      <- a SPIN-1 plan
tensor-flip        built OK, categories j = [1.0, 1.0]
helicity-flip      built OK, categories j = [1.5, 1.5]
```

and running the advised plan throws three frames down:
`RuntimeError: InclusiveKernel::amplitudes: spin state J = 1.000000 is not the
kernel's ion spin 1.500000`. The message names a plan that cannot work.

**F3 — `spin32_populations` refuses ordinary-looking (P_z, T).**
`lg.make_plan("helicity-flip", j=1.5, pz=0.7, pzz=0.6)` raises
*"unphysical (pz, t, o): negative population"*. Correct physics (the J = 3/2
Bloch domain is smaller than the spin-1 one at R₃ = 0), but the message does not
say which pair is reachable.

**F4 — a three-line Python segfault.** `ClusterPartialWave`'s `k` setter
(`python/bindings.cpp:1376-1378` (as of adec442)) calls `rebuild()`, which builds
`CubicSpline(k, phi)` (`src/core/b1_nuclear.cpp:177-180` (as of adec442)) with `phi` still
empty; the constructor indexes `y_[i+1]` (`src/core/b1_nuclear.cpp:66` (anchor read 2026-09-15 "d[i] = (y_[i")) out of
bounds. Reproduced:

```python
import lipolgen._lipolgen as _l
w = _l.ClusterPartialWave(); w.k = [0.1, 0.2, 0.3]     # Segmentation fault
```

Not a physics defect, but it is a crash reachable from the documented API, and
it is why §9's script goes through `from_uw` instead.

### 1.6 What is *not* broken

Everything downstream of the missing input works. With a caller-supplied kernel
carrying `b1_32_func = 0.05·F1` and `delta_32_func = −1e-2·F1`, on the same
hand-built ⁷Li A_zz plan, 40 000 events, seed 11:

```
sigma+ = 565389 pb   sigma- = 616516 pb   A_T = -0.043258
meta b1_model = 'caller-supplied kernel', b1_unpol = 'caller-supplied kernel'
```

−0.0433 against the −0.05/(1+εR) of `SPIN32_FINITE_GAMMA.md` Eq. (48) with
⟨εR⟩ ≈ 0.15. **The sampler, the tensor weight, the per-category cross sections,
the P8 spin guard and the provenance line are all already correct for J = 3/2.**
What is missing is a physics input and a first-class surface — not machinery.

---

## 2. Why ⁷Li's rank-2 input is *purely orbital*

### 2.1 The α + t channel is a single L = 1 wave — a selection rule, not a model

α is 0⁺ and the triton is ½⁺, so the pair's channel spin is S = ½ and the
parity of the relative motion is (−1)^L. ⁷Li's ground state is **3/2⁻**:

* parity ⇒ **L odd**;
* L ⊗ ½ ∋ 3/2 ⇒ **L ∈ {1, 2}**.

The intersection is **L = 1 alone**. L = 3 ⊗ ½ gives only 5/2 and 7/2, so there
is no D-wave/F-wave partner and **no interference term anywhere in the α–t
channel**. The repository already records this
(`docs/CONVENTIONS.md:249` (as of adec442), `docs/PHYSICS_CHANNELS.md:221` (as of adec442),
`docs/open_items/vmc_reconciliation.md:107-113`,
`tests/test_tagged.cpp:1198` — *"Aat11 is a selection-rule zero"*), and the data
files agree: `li7.at`'s second column `Aat11(k)` is the ½⁻ **excited** state's
amplitude and is MC noise at the 1e−4 level in the 3/2⁻ block;
`momenta/li7_at1.momentum` is the ½⁻ state's own full-size distribution
(S = 0.98683) and `momenta/li7_at3.momentum` is the 3/2⁻ ground state
(S_αt = **1.0084**, `VMC_S_ALPHA_T_LI7`, `include/lipolgen/tagged.hpp:109`).

**So the "L split" of the momenta block is a J split, not an L split.** There is
one wave.

### 2.2 Neither cluster carries rank-2

* α: J = 0 ⇒ no rank-1, no rank-2. (⁶Li's term (2α) exists for exactly the same
  reason and by the same argument: b₁^α ≡ 0, but its light-cone density carries
  the orbital alignment.)
* triton: J = ½ ⇒ **no rank-2 structure function exists for it.** Its internal
  structure enters b₁(⁷Li) only through the **spin-blind** F₁_t.

Therefore, in the two-cluster picture,

> **b₁(⁷Li) is 100 % orbital.** ⁶Li's term (1) — the embedded deuteron's own
> b₁_d, which is what `LI6_B1_RANK2_TRANSFER` = 0.921947 × (2/6) encodes — has
> **no ⁷Li counterpart**. Neither does term (3), the Clebsch–Gordan
> depolarisation of a spin-1 constituent, because there is no spin-1
> constituent.

### 2.3 P_zz for J = 3/2 in *this* repository

`spin.hpp:89-91` and `src/core/xsec.cpp:167-172`:

        J = 1    : P_zz = ⟨3J_z² − 2⟩ ∈ [−2, +1] ,   Q_NN(m) = (3m² − 2)/3
        J = 3/2  : T    = ⟨3J_z² − J(J+1)⟩/3 ∈ [−1, +1] ,  Q_NN(m) = (3m² − 15/4)/3 = m² − 5/4

so for J = 3/2 the alignment weight is **Q_NN(±3/2) = +1 and Q_NN(±½) = −1**,
and for a pure state Q_NN = T. `tensor_moments` returns the pair (Q_NN, 3·Q_NN)
for every spin: the b-sector multiplies Q_NN, the Δ (cos 2φ) sector multiplies
3·Q_NN. **`helicity_flip_plan`'s `pzz` field carries T, not P_zz, when j = 3/2**
(`include/lipolgen/bookkeeping.hpp:155` (anchor read 2026-09-15 "the normalized T in")); the third moment R₃ is set to
zero by `spin32_populations`'s default and there is no polarimeter for it
(`SPIN32_FINITE_GAMMA.md` §2.4). None of that reaches an unpolarised-beam
observable — the rank-≤2 theorem of `SPIN32_FINITE_GAMMA.md` §4 — so **A_zz and
cos 2φ on ⁷Li depend on T alone.**

The sign convention the slots must be filled in is fixed and is **not** the one
of arXiv:2209.12161: `SPIN32_FINITE_GAMMA.md` Eq. (23)–(24),

        F1^{(m)}|code = F1 − Q_NN(m) · b1_32|code ,   b1_32|code = − b1^{[Fu22]}

i.e. **b₁ > 0 means the |m| = 3/2 stretched doublet has the *smaller* cross
section.** Everything below is in the code's sign.

### 2.4 The rank-2 alignment of an L = 1 cluster wave — exact

Write the α–t relative-momentum density in ion state M. Only |φ₁(k)|² and the
Clebsch–Gordan weights enter, because there is one wave:

        n_M(k, k̂) = |φ₁(k)|² · A_M(k̂) ,
        A_M(k̂)   = Σ_{m_L m_S} |⟨1 m_L ; ½ m_S | 3/2 M⟩|² |Y_{1 m_L}(k̂)|²      (1)

Using ⟨L m|P₂(cos θ)|L m⟩ = [L(L+1) − 3m²]/[(2L−1)(2L+3)], which for L = 1 is
(2 − 3m_L²)/5, the sum collapses to a single Legendre coefficient:

        A_M(k̂) = (1/4π) [ 1 − Q_NN(M) · P₂(cos θ_k) ] ,   exactly, for every k   (2)

**Measured** (`clebsch_gordan` from the library, §9 block 6):

| M | ⟨P₂(cos θ_k)⟩ | −5⟨P₂⟩/Q_NN |
|---|---|---|
| ±3/2 | −0.200000 | **+1.000000** |
| ±½ | +0.200000 | **+1.000000** |

Two things follow at once.

**(a) Eq. (2) *is* the repository's ⁷Li polarimeter.** ⟨P₂⟩_M = −Q_NN(M)/5, and
for any fill ⟨P₂⟩ = −T/5 — which is precisely
`tests/test_tagged.cpp:383-400`, passing at tolerance 3e−4 today. The alignment
input to a ⁷Li b₁ is therefore **already gated by a committed test**, before any
b₁ exists.

**(b) The coefficient is k-independent.** A single partial wave cannot produce a
k-dependent alignment; ⁶Li's can and does, because its alignment is an S–D
*interference* whose sign flips at the S node (0.678 fm⁻¹) and again at the D
node (2.25 fm⁻¹) — `b1_nuclear.hpp:193-208`. **⁷Li has no node problem.**

### 2.5 The same derivation reproduces the code's own ⁶Li coefficients

Eq. (2) is not a new convention; it is the one `b1_nuclear.cpp:355-356` already
uses, transposed. Running the identical CG algebra on ⁶Li's α–d channel
(L = 2, S = 1, J = 1):

| channel | ⟨P₂⟩ at M | −5⟨P₂⟩/Q_NN |
|---|---|---|
| pure L = 2, M = ±1 | −0.100000 | **+1.500000** |
| pure L = 2, M = 0 | +0.200000 | **+1.500000** |
| pure L = 0 (any M) | 0 | — (an S wave is isotropic) |

and `b1_nuclear.cpp:356` writes `w_dd = 3·k·(W²/4)·(3c*²−1)`, i.e. a coefficient
of **3·2/4 = 1.500000** on `[W² P₂]`. Exact agreement.

The S–D term closes the loop too. Building the full angular density
A_M(k̂) = Σ_{m_S} |Σ_L φ_L(k) Y_{L m_L}(k̂) ⟨L m_L 1 m_S|1 M⟩|² with the CDKS
phases φ₀ = U, φ₂ = −W at U = W = 1 and projecting onto P₂ by quadrature on
200 001 points gives, at M = +1,

        a₂ = −1.91424      against      −Q_NN(1)·[ 6/√2 + 3/2 ] = −(1/3)(5.742641) = −1.914214

— the code's `w_sd` coefficient `3·(UW/√2)·(3c*²−1)` = (6/√2)·[UW·P₂] and its
`w_dd` = 1.5·[W²·P₂], both recovered from the same construction to five digits.

> **This is the load-bearing validation of §3's normalisation.** The ⁷Li
> coefficient 1.000 is produced by the identical procedure that reproduces the
> shipped ⁶Li coefficients 4.242641 and 1.500000.

---

## 3. The convolution kernel for a J = 3/2 target

### 3.1 The master formula

Let the CDKS δ-function fix cos θ* = c*(k, y) (`ConvolutionKinematics::cos_star`,
`src/core/b1_nuclear.cpp:281-285` (as of ac22331)) and write the repository's own prefactor
`pre = y·m_struck/(2κ)` (`src/core/b1_nuclear.cpp:377`). Inserting Eq. (2) into
the light-cone reduction gives, with no approximation beyond the ones the ⁶Li
kernel already makes,

        f_M(y) = pre ∫_{k_range(y)} k dk |φ₁(k)|² [ 1 − Q_NN(M) P₂(c*(k,y)) ]
               = f₁(y) − Q_NN(M) · f₁^{P₂}(y)                                    (3)

        f₁(y)     ≡ pre ∫ k dk |φ₁|²             (the unpolarised α–t density)
        f₁^{P₂}(y) ≡ pre ∫ k dk |φ₁|² P₂(c*)      (the tensor density)            (4)

Convolving per nucleon and matching to F1^{(M)} = F1 − Q_NN b1 (§2.3):

> **b₁(⁷Li)(x, Q²) = (3/7) ∫ (dz/z) f₁^{P₂,t}(z) F₁^t(x/z, Q²)
>                  + (4/7) ∫ (dz/z) f₁^{P₂,α}(z) F₁^α(x/z, Q²)**      (5)

with the two kinematics {m_struck, m_recoil} = {M_t, M_α} and {M_α, M_t},
ε = 2.467 MeV (`LI7_ALPHA_TAG().separation_energy`), and the counting factors
**3/7 and 4/7** replacing ⁶Li's 2/6 and 4/6. **Two terms, not four.**

Eq. (4) is *literally* `LightConeDensities::f_d` and `::f_d_p2`
(`include/lipolgen/b1_nuclear.hpp:411`) with the L = 1 wave placed in the φ₂
slot — because those two members carry the same `pre` and differ only by the
P₂(c*) weight (`src/core/b1_nuclear.cpp:351, 357, 360, 364`). That is how §5's
numbers were produced, and §7 records why it is an *abuse* of the interface
rather than an implementation.

### 3.2 What is exact and what is an ansatz

**Exact, given the α + t two-body picture:**

1. L = 1 is the only allowed relative wave (§2.1) — a selection rule.
2. b₁^t ≡ 0 and b₁^α ≡ 0 (§2.2) — a spin statement.
3. The alignment coefficient of Eq. (2), and hence the "1" in front of
   f₁^{P₂} in Eq. (3) — Clebsch–Gordan, verified two ways (§2.4, §2.5).
4. ⟨P₂⟩ = −T/5 for any fill — the same algebra, and a passing test today.
5. ∫ b₁ dx = 0 over the full support: ∫ f₁^{P₂} dy = 0 identically (the P₂ weight
   at fixed y integrates the sphere), so ∫dx ∫(dz/z) f^{P₂}(z) F₁(x/z) =
   [∫f^{P₂}dz][∫F₁du] = 0 **for any F₁**. Measured residuals on the shipped grid:
   **+5.0e−08 (struck t) and +4.7e−09 (struck α) against max|f^{P₂}| = 5.1 and
   6.8** — i.e. 1e−8 relative. See §6.4 for why this is a quadrature check and
   not a physics one.

**Ansatz, inherited from the ⁶Li design and *not* re-derived here:**

6. The plane-wave impulse/light-cone reduction itself: CDKS Eqs. (16), (17),
   (21) applied one cluster level up, with the "fair-share" variable
   z = (E − k_z κ)/m_struck (`b1_nuclear.hpp:294-329`, design 1.6). ⁷Li's mass
   ratios are milder than ⁶Li's (M_α/M_t = 1.327 against M_α/M_d = 1.987), so
   the documented 0.1–2.1 % cost of the mass-normalised z is, if anything,
   smaller — **not measured here**.
7. κ = 1 (CDKS Eq. (17) as printed) rather than the exact Eq. (21)
   √(1+γ²). Measured cost, §5.4: **+0.35 % to +69 %**, much larger than ⁶Li's
   −1 %…+7 %, because ⁷Li has no term (1) to dilute it.
8. Renormalising the light-cone density to S_αt = 1.0084. See §7 D3.
9. F₁ per nucleon of the α: isoscalar, no α EMC effect — carried over verbatim
   from ⁶Li (`Li6ConvolutionOptions::alpha_f1`, design A9).
10. b₂ = 2x·b₁ from the `TensorSF` base class, and Δ_32 = whatever is chosen.

**Wrong if carried over unexamined:**

11. **F₁ of the triton is not isoscalar.** ⁶Li's code sets `f1_alpha_ = f1_d_`
    and that is exactly right for a Z = N = 1 deuteron. The triton is Z = 1,
    N = 2. Measured effect of using `NuclearF2(TRITON())/3` in place of the
    isoscalar (F₁p+F₁n)/2 in the struck-t term: **+0.60 / +2.13 / +5.40 /
    +8.57 / +5.79 / −1.74 %** at x = 0.05 … 0.70.

### 3.3 What transfers from `b1_nuclear.hpp` and what does not

| ⁶Li machinery | Transfers to ⁷Li? |
|---|---|
| `ClusterPartialWave` | **No, as typed.** `from_vmc` and `from_uw` throw *"i^L is real for even L only"* (`src/core/b1_nuclear.cpp:210` (as of ac22331), `:233` (as of ac22331)). For a single wave the global phase is unobservable, so the refusal protects nothing here — but the type must be taught to say so |
| `ConvolutionKinematics` (`k_range`, `cos_star`, `y_max`) | **Yes, unchanged.** Only the two masses and ε change. Measured y_max = 1.6626 (struck t) and 1.3761 (struck α) |
| `LightConeDensities` quadrature (3-segment y grid, Simpson in k, exact k endpoints, baryon-number renormalisation) | **Yes, unchanged.** Converged: doubling every grid moves b₁ by ≤ 0.05 % |
| `LightConeDensities`' φ₀/φ₂ **interface** and its `delta_t_f` SD/DD split | **No.** There is no S–D interference to split. `f_d`/`f_d_p2` carry the needed objects but under names that lie about what they hold |
| `convolve(dens, g, x, x_max_g)` | **Yes, unchanged**, with `x_max_g = 1` because there is no injected b₁ table to run past x = 1 |
| Four-term structure (1)+(2d)+(2α)+(3) | **No — two terms.** (1) and (3) do not exist (§2.2) |
| `LI6_B1_PER_NUCLEON` = 2/6 and the ×2 on the α term | **Replaced** by 3/7 and 4/7 |
| `LI6_B1_RANK2_TRANSFER` = 0.921947 | **Does not transfer, and must not be reused.** It is 1 − (9/10)·P_D(α–d), the tensor dilution of an *embedded spin-1 deuteron* inside an L = 2 component. ⁷Li has no embedded spin-1 cluster. The number that plays the analogous role is the −1/5 of ⟨P₂⟩ = −T/5, and it is exact rather than model-dependent |
| Phase A's **normalisation convention** (per-nucleon vs per-deuteron, `B1_CDKS_TABLE_TO_PER_NUCLEON` = 1, `B1_MILLER_TABLE_TO_PER_NUCLEON` = 0.5) | **Does not arise.** Those constants exist to reconcile two published *deuteron b₁ tables*. There is no ⁷Li table and no A = 3 table to reconcile. What replaces it is the counting choice 3/7 + 4/7 and the F₁ normalisation of §3.2 item 11 |
| Phase A's `r_func` decision (`r_sigma_lt`, not `r1998`, `b1_nuclear.hpp:487-507`) | **Transfers, and matters more.** The reason given there — the orbital terms convolve F₁ against a density that integrates to zero and so respond to the *slope* of R — applies to **all** of b₁(⁷Li). Measured: −9.66 % … −46.02 % (§5.4) |
| The **mandatory ±100 % band** (`banded()`, from Q(⁶Li) vs Q_d) | **The justification does not transfer** — §6.2 shows ⁷Li's quadrupole is *predicted to 13 %*, not cancelled. But §5.4 shows a **larger** band is needed for a different reason: the unpolarised input |

---

## 4. What can be gated at A = 7 — and why there is no A = 3 gate

### 4.1 There is no A = 3 analogue of the A = 2 gate, and the reason is a selection rule

The A = 2 gate works because the same convolution kernel, run one cluster level
down (p + n instead of α + d), produces b₁ of a nucleus for which a published
theory curve exists (CDKS Fig. 4). The A = 3 analogue of α + t would be the
triton as t = d + n or ³He = d + p — and **³H and ³He are J = ½.** A spin-½
target has no rank-2 multipole and therefore **no b₁ at all**, published or
otherwise. `TaggedModel::tensor_dilution` throws for exactly this reason
(`src/core/tagged.cpp:527-530`, *"tensor dilution defined for S_c = 1"*).

**So the escalation clause of design §5.4 has no A = 3 form. Nothing about the
⁷Li convolution can be validated by running it on a lighter system.** That must
be said in any header that ships this model, and it is the single sharpest
difference from the ⁶Li situation.

### 4.2 The gate that *does* exist: the ⁷Li quadrupole moment

Both α–t clusters have zero intrinsic quadrupole (0⁺ and ½⁺), so — exactly as
for b₁ — **the whole of Q(⁷Li) is orbital in this picture**, and it is the same
rank-2 alignment. For a pure L = 1, S = ½, J = 3/2 two-cluster state the
spectroscopic quadrupole factorises:

        Q = Z_eff ⟨r²⟩ ⟨3cos²θ − 1⟩_{M=3/2} = − (2/5) Z_eff ⟨r²⟩                (6)
        Z_eff = Z_α (M_t/M₇)² + Z_t (M_α/M₇)² = 0.695075   (A-number form 34/49 = 0.693878)

with ⟨3cos²θ−1⟩_{M=3/2} = 2[L(L+1) − 3m_L²]/[(2L−1)(2L+3)] = −2/5 at m_L = +1 —
the *same* −1/5 that Eq. (2) puts in front of P₂, times 2.

Computed from `momenta/li7_at3.momentum` by Fourier–Bessel
(u₁(r) = r√(2/π)∫k²j₁(kr)φ₁(k)dk; the transform round-trips the k-space norm to
3e−6):

```
<r^2> = 12.5348 fm^2   (r_rms = 3.5405 fm)
Q(7Li)_alpha-t = -3.4851 fm^2
   grid/truncation: -3.4720 (r_max 20) / -3.4851 (30) / -3.4858 (40) / -3.4862 (60)
   VMC MC band:     -3.4964 (-1 sigma) / -3.4851 (0) / -3.4737 (+1 sigma)
```

Against the measured moment — which **was** an external number absent from
this tree when this note was written (§7 D11), and **since 2026-09-06 is
`LI7_QUADRUPOLE_FM2` = −4.06 fm² in `rc.hpp`** — Q(⁷Li):

> **ratio 0.8584** against the committed −4.06 fm² (TUNL A = 5, 6, 7, the same
> evaluation `LI6_QUADRUPOLE_FM2` comes from), **0.8713** against the −4.00(3)
> fm² this note originally quoted (Voelk *et al.*, NPA 530 (1991) 475, one of
> the eight determinations N. J. Stone lists without recommending one).
> **13–14 % low, converged to 0.1 %, MC band ±0.3 %.**

**This is no longer an offline calculation.** `li7_alpha_t_quadrupole`
(`b1_nuclear.hpp`) ships it, and it reproduces §9's `q_at` **to the last bit**
(relative difference 0.000e+00 on ⟨r²⟩, Q and Z_eff, measured 2026-09-06); the
grid row, the MC band and both ratios are pinned by `tests/test_b1_nuclear.cpp`
T13 and `python/tests/test_li7_rank2.py` G8, and the run-2026-09-06 numbers are
in `../run_2026-09-06/phase_B_numbers.md` §B3. What the test asserts is the
**ratio**, reported — never a pass/fail on b₁(⁷Li), which is still not
implemented (open item 15) and which G8 re-checks is still exactly zero in the
very test that runs the gate.

For contrast, from the same code path on ⁶Li:
`alpha_d_quadrupole_fm2(li6_alpha_d_partial_waves())` = **−0.3333 fm²** against
`LI6_QUADRUPOLE_FM2` = **−0.0818 fm²** — a **factor 4.08**, which is why the ⁶Li
design could only ever claim it as a *sign* gate.

**This is the A = 7 gate, and it is a genuine one:** a measured rank-2
observable of the same nucleus, predicted with no free parameter by the same
wave function and the same CG coefficient the b₁ construction uses, to 13 %.
It is offline-doable today (§9 reproduces it in ~20 lines) and it should be
committed **whether or not the b₁ is**, because it validates the *input* that
`li7_vmc_waves` already ships to the tagged channel.

Its honest limits: (i) the α–t overlap is not the full ⁷Li wave function —
S_αt = 1.0084 > 1 is itself proof that it is not a probability — and the missing
13 % is the expected size of cluster distortion and non-α–t components; (ii) a
13 % agreement gates the wave function's ⟨r²⟩ and its P-wave character, **not**
the light-cone convolution or the DIS input, which is where §5.4's real
uncertainty lives.

### 4.3 The ⟨P₂⟩ = −T/5 polarimeter

Already a passing committed test (`tests/test_tagged.cpp:336-353`, tolerance
3e−4, both pure states and two mixed fills). §2.4 shows it is the *same*
Clebsch–Gordan coefficient that Eq. (3) needs. So the alignment half of the
construction is gated today.

What it does **not** gate: the radial wave function (it is a ratio, so |φ₁|²
cancels), the light-cone reduction, the counting factors, or F₁. It is
necessary, not sufficient.

### 4.4 Sum rules and positivity — what they test

* **Close–Kumano ∫b₁ dx = 0.** For a two-cluster convolution it is *automatic*
  (§3.2 item 5) and holds for any F₁, so it tests the **quadrature and the y
  support**, not the model. This is the same status the ⁶Li code gives it
  (`close_kumano_integral` — "REPORTED, never enforced"). On the truncated
  window [0.01, 1.2] the ⁷Li estimate gives ∫b₁ dx = **+7.51e−05** against
  max|b₁| = 1.11e−03, and ⁶Li's own reports **+1.459e−04** — both are the
  truncation, not the physics.
* **[Fu22] Eq. (21) ∫g₁_rank3 dx = 0.** Concerns the rank-3 slot the code
  deliberately does not have (`xsec.hpp:220-224`). Not a b₁ gate.
* **J = 3/2 positivity, `SPIN32_FINITE_GAMMA.md` Eq. (27):**
  |(m/J)g₁ + R₃(m)g₁_rank3| ≤ F₁ − Q_NN(m)·b₁_32. With the measured
  b₁/F₁ ≤ 1.04e−02 this is satisfied by five orders of magnitude and is not a
  discriminating test at these sizes. It *would* become one for any toy input of
  order 0.05·F₁ — which is exactly the reference scenario
  `validation/dump_polligen_reference.py:491` (anchor read 2026-09-15 "0.05 * f1") uses, so it is worth having, but
  it does not gate a physical b₁(⁷Li).
* **The ⁶Li-core cross-check** (§5.5) is not a gate — it is a second model.

**Summary: the only *external, measured, discriminating* offline gate available
for ⁷Li's rank-2 input is the quadrupole moment.**

---

## 5. The numbers

All at Q² = 2.5 GeV², per nucleon, in the code's sign convention (§2.3).
Nominal configuration: `ToyF2` (the shipped default), `r_sigma_lt`, CDKS
Eq. (17) with κ = 1, S_αt = 1.0084, F₁ per nucleon isoscalar for both clusters,
quadrature = `LightConeDensities::Options()` defaults.

### 5.1 The leading estimate

```
x                                          0.05         0.10         0.20         0.30         0.50         0.70
b1(7Li) per nucleon                  1.2463e-04   5.7001e-05   3.7480e-05   3.3103e-05   5.4849e-05   1.0400e-04
x b1                                 6.2313e-06   5.7001e-06   7.4959e-06   9.9308e-06   2.7424e-05   7.2797e-05
b1/F1  (the kernel's own F1)          3.576e-05    3.848e-05    6.899e-05    1.294e-04    8.944e-04    1.044e-02
```

Split between the two terms: struck-t : struck-α = **1 : 0.756**, against the
analytic counting-and-mass expectation [(4/7)/(3/7)]·(M_t/M_α)² = **0.757** —
the ⁷Li form of the consistency check `b1_nuclear.hpp` states for ⁶Li as
2(M_d/M_α)² = 0.5064.

The observable, from `SPIN32_FINITE_GAMMA.md` Eq. (48):

        A_T^{(3/2)}(θ_S = 0) · (1 + εR) = − b1_32/F1

so the alignment contrast between the T = +1 and T = −1 fills is **−3.6e−05 at
x = 0.05, rising to −1.0e−02 at x = 0.70**, before the 1/(1+εR) factor. That
factor is kinematics-dependent and is *not* computed per x here; the only value
measured in this note is the sample average of the §1.6 run, 0.043258/0.05 =
**0.865** at the default ⁷Li inclusive configuration. Note the
**2/3 is absent** — Eq. (48)'s warning: quoting −(2/3)b₁/F₁ for a contrast that
returns −b₁/F₁ is a 50 % error. §7 D6.

### 5.2 Against ⁶Li

```
x                                          0.05         0.10         0.20         0.30         0.50         0.70
7Li  alpha-t (this note)             1.2463e-04   5.7001e-05   3.7480e-05   3.3103e-05   5.4849e-05   1.0400e-04
6Li  orbital (2d)+(2a)               4.1581e-05   1.9053e-05   1.2529e-05   1.1075e-05   1.8411e-05   3.4763e-05
6Li  convolution total, 4 terms      7.3014e-05  -2.7966e-05  -1.0603e-04  -1.3911e-04   1.0348e-04   4.1097e-04
6Li  SHIPPED Li6B1(MillerB1)         7.3006e-03   4.1222e-03   1.3350e-03   2.2389e-05  -3.7300e-04   1.7232e-04
```

Three statements, each with its window:

1. **⁷Li's b₁ is 2.99 × ⁶Li's orbital term**, at every x tested and under every
   input swap:

   | configuration | 0.05 | 0.10 | 0.20 | 0.30 | 0.50 | 0.70 |
   |---|---|---|---|---|---|---|
   | ToyF2 | 2.9972 | 2.9917 | 2.9915 | 2.9889 | 2.9792 | 2.9916 |
   | MSTW2008 LO | 3.0160 | 2.8623 | 2.9772 | 3.0013 | **5.4390** | 2.9950 |
   | ToyF2 + r1998 | 2.9979 | 2.9927 | 2.9889 | 2.9874 | 2.9760 | 3.0029 |
   | ToyF2, finite-\|q\| δ | 2.9931 | 2.9920 | 2.9917 | 2.9881 | 2.9726 | 2.9901 |

   **The single exception is the boxed 5.44**, where ⁶Li's orbital term under
   MSTW is −5.29e−07, i.e. crossing zero; the ratio there is a division by
   nearly nothing, not a physics statement. Excluding it, the spread is
   **2.86 … 3.02**.

2. **Against ⁶Li's convolution total**, ⁷Li is comparable in magnitude
   (factor 0.24 … 1.7) and **single-signed positive over 0.05 ≤ x ≤ 0.70**,
   where ⁶Li's changes sign twice. The sign difference is not a ⁷Li statement:
   it is ⁶Li's term (1) (the CDKS b₁_d camp) dominating and changing sign.

3. **Against what a ⁶Li user actually gets today** — `Li6B1(MillerB1)`, the
   shipped default — ⁷Li's estimate is **1.71 % / 1.38 % / 2.81 %** of it at
   x = 0.05 / 0.10 / 0.20. Above x = 0.20 the ratio stops meaning anything,
   because Miller's own curve crosses zero near x = 0.30 (+2.24e−05 there) and
   is negative at x = 0.50: the ratios are +148 %, −15 %, +60 % at
   x = 0.30 / 0.50 / 0.70. That gap is the Miller-versus-CDKS **camp** gap
   (`PHYSICS_CHANNELS.md:184`: two orders of magnitude and a sign), not a
   ⁶Li/⁷Li physics gap, and it must never be quoted as one.

### 5.3 Why 2.99, and why it is F₁-independent

The tensor kernel integrates to zero, so ∫(dz/z)f(z)F₁(x/z) has no
leading term; expanding F₁(x/z)/z about z = 1, the convolution is controlled by
the kernel's **first non-vanishing moment** M_n = ∫(z−1)^n f(z)dz. Measured,
with the counting factors already applied:

```
7Li t  counting 0.4286  M0 +4.999e-08  M1 +7.835e-04  M2 +7.698e-04
7Li a  counting 0.5714  M0 +4.731e-09  M1 +4.448e-04  M2 +4.373e-04
6Li d  counting 0.3333  M0 +3.411e-08  M1 +3.924e-04  M2 +3.857e-04
6Li a  counting 0.6667  M0 -3.660e-10  M1 +9.932e-05  M2 +9.775e-05
weighted:  M0  7Li +2.4127e-08  6Li-orb +1.1127e-08   (both zero to quadrature)
           M1  7Li +5.8995e-04  6Li-orb +1.9701e-04   ratio 2.9945
           M2  7Li +5.7981e-04  6Li-orb +1.9374e-04   ratio 2.9927
```

**M₁'s ratio is 2.9945 — the observed b₁ ratio.** The two kernels have quite
different y-shapes (checked: their pointwise ratio swings from +0.09 to −174
across y ∈ [0.8, 1.2] and changes sign), so this is *not* a shape identity; it
is the leading moment coinciding. That is also the whole explanation of the
F₁-independence: any F₁ sees the same first moment.

The naive k-space alignment ratio is **not** 3 and must not be quoted:
∫k²[(6/√2)UW + 1.5W²]dk / ∫k²(U²+W²)dk = **0.120433** for ⁶Li (SD 0.0914,
DD 0.0290) against **1.000000** for ⁷Li — a factor 8.3, which the P₂(cos θ*)
correlation and the different y widths (⟨k⟩ = 0.1860 GeV for α–t against
0.1224 GeV for α–d) cut to 3.

### 5.4 Systematics — and the sign is *not* predicted

Percentage change of b₁(⁷Li) against the nominal row:

```
x                                          0.05         0.10         0.20         0.30         0.50         0.70
R -> r1998                               -9.66%      -36.10%      -46.02%      -33.48%      -10.76%       +4.74%
finite-|q| delta (CDKS Eq. 21)           +0.35%       +1.42%       +5.67%      +12.78%      +36.07%      +69.32%
S_at -> 1                                -0.83%       -0.83%       -0.83%       -0.83%       -0.83%       -0.83%
unpol -> MSTW2008 LO                    -79.35%     -123.38%     -237.77%     -347.64%     -105.24%      +79.10%
VMC MC -1 sigma                          -0.42%       -0.42%       -0.42%       -0.43%       -0.44%       -0.43%
VMC MC +1 sigma                          +0.45%       +0.45%       +0.46%       +0.46%       +0.48%       +0.46%
quadrature x2                            -0.05%       -0.05%       -0.02%       -0.01%       -0.00%       -0.00%
true F1_t (Z=1,N=2) not isoscalar        +0.60%       +2.13%       +5.40%       +8.57%       +5.79%       -1.74%
```

Read this table in the order the numbers are small to large.

* **Numerics are not the problem.** Doubling every grid: ≤ 0.05 %. The VMC file's
  own MC error band: ±0.5 %. The spectroscopic-factor convention: 0.83 %.
* **The model choices the ⁶Li design already agonised over cost tens of percent
  here** — R (up to 46 %), κ (up to 69 %), the triton's isospin (up to 8.6 %).
  Each is larger than its ⁶Li counterpart for the same reason: ⁶Li's term (1)
  carries an injected b₁_d table with no F₁ and no R in it, so it dilutes them.
  **⁷Li has nothing to dilute them with.**
* **The unpolarised backend flips the sign.** `ToyF2 → MSTW2008 LO` gives
  b₁(⁷Li) = +2.573e−05, **−1.333e−05, −5.164e−05, −8.198e−05, −2.877e−06**,
  +1.863e−04 at x = 0.05…0.70: **negative at four of six points where the
  nominal row is positive at all six.**

  This is not a ⁷Li pathology. The identical swap on ⁶Li's *orbital terms alone*
  gives MSTW/Toy = **0.205, −0.244, −1.384, −2.466, −0.029, +1.789** — the same
  behaviour, hidden in the shipped ⁶Li answer by term (1). The A = 2 gate's own
  verdict says MSTW is the configuration to quote (`b1_nuclear.hpp:25-40`,
  condition 1). **So the configuration the repository tells you to publish is
  the one that reverses the sign of the ⁷Li prediction.**

> **The measured claim: with the inputs committed to this tree, an α–t
> convolution predicts |b₁(⁷Li)/F₁| of order 3e−05 … 1e−02 across
> 0.05 ≤ x ≤ 0.70 at Q² = 2.5, and does *not* predict its sign.**

### 5.5 A second in-tree decomposition disagrees — order of magnitude *and* sign

`data/vmc/li7_li6_n/li6n_31.table` is the VMC ⁷Li → ⁶Li(1⁺) + n overlap,
p3/2 and p1/2 channels with S = 0.436(2) and 0.246(1). In *that* picture ⁷Li's
core is **spin 1 and has a b₁ of its own** — the analogue of ⁶Li's term (1) that
α + t lacks. The Clebsch–Gordan transfer (this is a **sketch**: the CG factor
only, with the orbital smearing and the struck neutron omitted, and the table is
r-space with no k-space partner in the tree):

```
channel j=3/2  S=0.436  <P_zz(6Li core)> at M=3/2 = -0.800   transfer (6/7)<P_zz>/3 = -0.228571
channel j=1/2  S=0.246  <P_zz(6Li core)> at M=3/2 = +1.000   transfer                = +0.285714
S-weighted b1(7Li)/b1(6Li) per nucleon = -0.043067        (total S = 0.682)
```

Note the **82 % cancellation** between the two neutron channels — a second
instance of the same near-cancellation that makes ⁶Li's quadrupole small.
Folding it onto the two ⁶Li b₁ camps:

```
  x    6Li-core route b1/F1: Miller camp   CDKS camp    |  alpha-t route (this note)
 0.05        -8.9971e-05   -8.9981e-07    |   +3.5760e-05
 0.10        -1.1919e-04   +8.0860e-07    |   +3.8483e-05
 0.20        -1.0461e-04   +8.3078e-06    |   +6.8990e-05
 0.30        -3.7007e-06   +2.2993e-05    |   +1.2939e-04
 0.50        +2.5333e-04   -7.0279e-05    |   +8.9445e-04
 0.70        -7.0685e-04   -1.6858e-03    |   +1.0436e-02
```

**Three estimates of the same quantity, spanning two orders of magnitude and
both signs at every x below 0.5.** They are not additive — S_αt = 1.008 and
S(⁶Li+n) = 0.682 are non-orthogonal decompositions of one nucleus, and summing
them would double count. §7 D1.

---

## 6. Size: is b₁(⁷Li) comparable to ⁶Li's, or suppressed?

**Neither "comparable" nor "suppressed" is the right word without a referent.
Measured, with the window:**

* **Against ⁶Li's *orbital* term — the like-for-like comparison, because that is
  all ⁷Li has — it is 2.99 ± 0.02 times LARGER** over 0.05 ≤ x ≤ 0.70 and under
  every input variation tried (one exception, §5.2). The intuition that "the
  triton carries no tensor, so only the P-wave alignment survives, so it must be
  small" is **wrong**: ⁶Li's orbital alignment is bought with a small D-state
  amplitude *and* partially cancelled across the S node; ⁷Li's is the whole wave
  function, coherent, with no node. The k-space alignment strength is 1.000
  against ⁶Li's 0.120 — a factor 8.3 (cut to 3 by the light-cone geometry, §5.3).
* **Against ⁶Li's convolution total**, 0.24 … 1.7 — the same order.
* **Against what ⁶Li users get today** (`Li6B1(MillerB1)`), 1.4 % … 2.8 % over
  x ≤ 0.20, and undefined above it because Miller's own curve changes sign — but
  that ratio is the Miller-vs-CDKS **camp** gap, not a ⁷Li statement.
* **In observable terms**, |A_T| ≈ 0.004 % at x = 0.05 rising to ≈ 1.0 % at
  x = 0.70, before 1/(1+εR).
* **The ⁶Li-core decomposition disagrees**, giving |b₁/F₁| up to 8× the α–t
  value with the opposite sign (§5.5).

So: **the correct answer to "how big" is that the leading estimate is of the
same order as ⁶Li's convolution, that its most likely honest presentation is a
band containing zero and both signs, and that the model spread — not the
statistics, not the quadrature, not the VMC errors — is what sets it.**

---

## 7. Author decisions this raises

| # | Decision | What is at stake |
|---|---|---|
| **D1** | **Which decomposition.** α + t (S = 1.008, purely orbital, this note) or ⁶Li + n (S = 0.682, carries a core b₁), or a refusal to choose | §5.5: they disagree in sign and by up to 8×. They are non-orthogonal; summing them double counts. α + t is the one already in the C++ (`li7_vmc_waves`, `li7_alpha_channel`) and the one with a k-space table and a quadrupole gate |
| **D2** | **The unpolarised backend the convolution folds against.** `ToyF2` (shipped ⁶Li default) or MSTW2008 LO (what the A = 2 gate says to quote) | **Decides the sign of b₁(⁷Li)** at four of six x points (§5.4). This is the same unresolved question as ⁶Li's `--b1-unpol`, but with no term (1) to hide it. It should be settled once, for both isotopes, not twice |
| **D3** | **The normalisation target.** S_αt = 1.0084 > 1 means the α–t overlap is **not a probability**, so the ⁶Li rule ("the non-α–d 18 % of ⁶Li is given b₁ = 0", `use_spectroscopic_factor`) has no ⁷Li form: renormalising to S puts *more* than all of ⁷Li in the α–t channel | Numerically 0.83 %, so this is a **provenance** decision, not a magnitude one — but it must be recorded, because "renormalise to the file's own S" and "renormalise to 1" are different claims about ⁷Li |
| **D4** | **F₁ of the triton.** Isoscalar (the ⁶Li shortcut) or the true Z = 1, N = 2 triton | 0.6 … 8.6 % (§3.2 item 11). The ⁶Li code's `f1_alpha_ = f1_d_` is *correct* for a deuteron and *wrong* here; copying it would be a silent error, not a documented approximation. `include/lipolgen/triton_sf.hpp` already exists for the tagged channel |
| **D5** | **κ: CDKS Eq. (17) or Eq. (21).** ⁶Li defaults to Eq. (17) because the switch costs −1…+7 % there | Costs +0.35 … +69 % here (§5.4). The ⁶Li default's stated justification ("the term κ multiplies is the small ORBITAL one") **inverts** for ⁷Li, where it is the only term |
| **D6** | **Which A_zz.** `A_T` (= −b₁/F₁) or `A_zz^{(3/2)} ≡ (2/(3T))(σ_T/σ_U − 1)` (= −(2/3)b₁/F₁) | A factor 3/2. `SPIN32_FINITE_GAMMA.md` §5.4 says either is defensible and §6.4 test 10 pins whichever is adopted, but the spin-1 `azz()`'s explicit `2.0/3.0` (`src/core/asymmetries.cpp:85` (anchor read 2026-09-15 "TENSOR_LL_SIGN")) has **no J = 3/2 counterpart**. Must be chosen before any ⁷Li tensor number is published |
| **D7** | **`ClusterPartialWave` refuses odd L** (`b1_nuclear.cpp:244`, `:267`) | For a single wave the i^L phase is unobservable, so the refusal protects nothing — but it must be *taught* that, in the type, rather than worked around by passing `l = 0` as §9 does |
| **D8** | **`LightConeDensities`' interface.** Its φ₀/φ₂ + SD/DD shape is an A = 2/⁶Li shape. §9 reuses `f_d`/`f_d_p2` with the L = 1 wave in the φ₂ slot — numerically exact, semantically a lie | A first-class implementation needs either an L-generic alignment slot or an explicit `Options::alignment_coefficient` (1 for L = 1/S = ½/J = 3/2, 1.5 for the L = 2 DD term, 6/√2 for the SD one) — which would also let the ⁶Li coefficients be *derived* rather than hard-coded at `b1_nuclear.cpp:355-356` |
| **D9** | **b₂_32.** `TensorSF`'s base gives 2x·b₁ | No reason to differ, but the ⁷Li slot is separate (`b2_32_func`) and silence there means 2x·b₁ by default, which should be stated rather than inherited |
| **D10** | **Δ_32 (cos 2φ gluon transversity).** ⁶Li gets `toy_delta_gluon(…, 1e-2)`; ⁷Li gets nothing | Filling `b1_32_func` alone leaves cos 2φ at zero. There is **no ⁷Li Δ model**; reusing the ⁶Li toy means adopting an arbitrary 1e−2 scale for a second nucleus. Note 3·Q_NN = ±3 for J = 3/2 against ±1/−2 for spin 1, so the same Δ gives a *larger* ⁷Li cos 2φ amplitude |
| **D11** | ~~**Q(⁷Li) is not in the tree.**~~ **PAID 2026-09-06**: `LI7_QUADRUPOLE_FM2` = **−4.06 fm²** now sits beside `LI6_QUADRUPOLE_FM2` at `rc.hpp:669` | Exactly one copy, from a real source — TUNL's A = 5, 6, 7 evaluation (NPA 708 (2002) 3), the **same document** the ⁶Li constant comes from, which is what fixed the choice; the −4.00(3) vs −4.06 spread is recorded in the constant's own comment together with all eight of Stone's ⁷Li entries, and the predicted price is confirmed: ratio **0.871265 → 0.858389** |
| **D12** | **The run-plan surface.** No spin-3/2 tensor plan exists (§1.5 F1) | A ⁷Li A_zz programme needs `tensor_thirds_plan`'s J = 3/2 analogue — and §2.3's note that for J = 3/2 the pure-alignment fill (P_z = 0, symmetric populations) is *rank-3 clean by construction*, so the honest ⁷Li A_zz plan is the two-state T = +1 / T = −1 contrast of §1.3, not a thirds pattern |

---

## 8. Recommendation

**Path (c): a three-step order, not (a) and not (b).**

### Step 1 — the loud zero and the four defects, now, independent of any physics

None of this needs a b₁ and none of it can be wrong.

1. **`meta["b1_model"]`** must not say `"miller"` on a run where no rank-2 slot
   was filled. Say what happened: `"none (spin 3/2: no rank-2 input)"`, and the
   same for `b1_band_scale` / `b1_alpha_d_dwave_weight` / `b1_unpol`. This is the
   repository's own rule (`src/core/pipeline.cpp:531-546` (as of adec442)) applied to the model
   name.
2. **A banner line on every ⁷Li inclusive run**, unconditional, in the shape the
   `--b1-model` block already uses: *"7Li rank-2: b1_32 = b2_32 = Delta_32 = 0 —
   the tensor and cos 2φ sector of this run is identically zero (no published
   b1 for 7Li; see docs/…/phase_D_li7_rank2.md)."*
3. **Fix `make_plan`'s advice** (§1.5 F2): at J = 3/2 only `helicity-flip`
   works, and it is not a tensor plan.
4. **Fix the `ClusterPartialWave` segfault** (§1.5 F4) — three lines of Python.

### Step 2 — commit the A = 7 quadrupole gate, before any b₁

`Q(⁷Li)_α–t = −3.4851 fm²` against the measured −4.00(3) fm² validates the
**wave function that `li7_vmc_waves` already ships to the tagged channel**, and
does so an order of magnitude better than the ⁶Li sign gate validates its own.
It is worth having on its own merits (it is the only quantitative check the
⁷Li α–t input has ever had), it costs one constant (D11) and ~20 lines, and it
is the thing §4.1's absence of an A = 3 gate makes indispensable.

**DONE 2026-09-06** (run 2026-09-06 task B3): the constant is
`LI7_QUADRUPOLE_FM2` = −4.06 fm² at `rc.hpp:669`, the ~20 lines are
`alpha_t_quadrupole` / `li7_alpha_t_quadrupole` in `b1_nuclear.{hpp,cpp}`, and
the gate is `tests/test_b1_nuclear.cpp` **T13** + `python/tests/test_li7_rank2.py`
**G8**. Step 2 is closed; steps 1, 3 and 4 of this section are not, and the
gate deliberately licenses none of them.

### Step 3 — the α–t b₁ as a third `B1Model`, but only after D2

Everything else is ready: Eq. (5) is 40 lines on top of `LightConeDensities`;
`b1_32_func` is already a live slot reachable from Python (§1.6); the sampler,
the P8 guard, the meta and the CLI plumbing all work at J = 3/2 today.

**What is not ready is the answer.** Shipping a b₁(⁷Li) whose *sign* depends on
whether the user typed `--b1-unpol toy` or `--b1-unpol mstw` — where the second
is the one the A = 2 gate tells them to use — would put the repository in the
position of publishing a tensor asymmetry whose direction is a flag. No band
expresses that: the band would have to contain zero and both signs, at which
point the number carries no information a loud zero does not.

So: **(a) is premature and (b) is insufficient.** Do (b) immediately because it
is a defect fix, do the gate because it is free and it is the only real
validation ⁷Li will get offline, and hold the b₁ until D2 is decided once for
both isotopes. When it is decided, this note's §5.1 is the number, §5.4 is its
band, §4.2 is its gate, and §7 D3–D10 are the eight lines of provenance that
have to ship beside it.

**One thing not to do:** reuse `LI6_B1_RANK2_TRANSFER` = 0.921947 for ⁷Li. It is
the tensor dilution of an *embedded spin-1 deuteron*; ⁷Li has no spin-1
constituent, and the number that plays its role — the −1/5 of ⟨P₂⟩ = −T/5 — is
exact rather than model-dependent (§3.3).

---

## 9. Reproduction

From a clean tree at `903fcc9`, `cd` to the repo root and `source env.sh`. The
script below regenerates every number in §§1.2, 2.4, 2.5, 3.2 item 5, 4.2, 5.1,
5.2, 5.3, 5.4 and 6. It writes nothing and mutates nothing.

```python
# scratch/li7_rank2_estimate.py  --  READ-ONLY
import math, numpy as np, lipolgen._lipolgen as _l

HB   = _l.HBARC_GEV_FM
UNIT = 1.0 / math.sqrt(2.0 * math.pi**2 * HB**3)          # b1_nuclear.cpp:257-258
MT, MA, M7 = _l.nuclear_mass(1,3), _l.nuclear_mass(2,4), _l.nuclear_mass(3,7)
MD   = _l.nuclear_mass(1,2)
EPS7 = _l.li7_alpha_channel().base.separation_energy      # 0.002467 GeV
EPS6 = _l.li6_alpha_channel().base.separation_energy
Q2, XS = 2.5, [0.05, 0.10, 0.20, 0.30, 0.50, 0.70]
v7 = _l.vmc_from_momentum(_l.data_path("vmc/momenta/li7_at3.momentum"), 0, 0, 1, None)

def waves(v):                       # from_uw refuses odd L; only |phi|^2 is used
    k = list(np.asarray(v.k, float)); p = list(np.asarray(v.psi, float) * UNIT)
    return (_l.ClusterPartialWave.from_uw(k, p, 0),
            _l.ClusterPartialWave.from_uw(k, [0.0]*len(k), 0))

def dens7(v, ms, mr, kappa=1.0, norm=_l.VMC_S_ALPHA_T_LI7, quad=None):
    w, z = waves(v)
    o = quad or _l.LightConeDensities.Options()
    o.renormalize = True; o.norm_target = norm
    k = _l.ConvolutionKinematics()
    k.m_struck, k.m_recoil, k.separation, k.kappa = ms, mr, EPS7, kappa
    return _l.LightConeDensities(z, w, k, o)              # L=1 wave in the phi2 slot

def f1_of(unpol, rf, ion=None, per=2.0):                  # CDKS Eq. (22), per nucleon
    f2 = _l.NuclearF2(ion or _l.deuteron(), unpol, None, rf)
    return lambda x: 0.0 if not (0.0 < x < 1.0) else _l.f1_cdks(f2.f2a(x,Q2)/per, x, Q2, rf)

def b1_li7(v=v7, unpol=None, rf=None, kappa=1.0, norm=_l.VMC_S_ALPHA_T_LI7, quad=None):
    f1, out = f1_of(unpol or _l.ToyF2(), rf), []
    for x in XS:
        kap = math.sqrt(1.0 + _l.gamma_squared(x, Q2)) if kappa == "fq" else kappa
        dt, da = dens7(v,MT,MA,kap,norm,quad), dens7(v,MA,MT,kap,norm,quad)
        out.append((3.0/7.0)*dt.convolve(dt.f_d_p2, f1, x, 1.0)      # Eq. (5)
                 + (4.0/7.0)*da.convolve(da.f_d_p2, f1, x, 1.0))
    return np.array(out)

# --- the quadrupole gate, Eq. (6) -----------------------------------------
def j1(x):
    x = np.asarray(x, float); o = np.empty_like(x); s = x < 1e-6
    o[~s] = np.sin(x[~s])/x[~s]**2 - np.cos(x[~s])/x[~s]; o[s] = x[s]/3.0; return o
def q_at(v, rmax=30.0, nr=4000, nk=8001):
    kgv = np.asarray(v.k, float)
    w = _l.ClusterPartialWave.from_uw(list(kgv), list(np.asarray(v.psi, float)), 0)
    kf = np.linspace(0.0, kgv[-1]/HB, nk); pf = np.array([w(k*HB) for k in kf])
    r  = np.linspace(rmax/nr, rmax, nr)
    u  = np.array([rr*math.sqrt(2/math.pi)*np.trapz(kf**2*j1(kf*rr)*pf, kf) for rr in r])
    n  = np.trapz(u*u, r); r2 = np.trapz(r*r*u*u, r)/n
    Z  = 2.0*(MT/M7)**2 + 1.0*(MA/M7)**2
    return r2, -0.4*Z*r2, Z

# --- the CG coefficient, and the check against the code's 6Li DD term ------
def p2LM(L, m): return (L*(L+1) - 3.0*m*m)/((2.0*L-1.0)*(2.0*L+3.0))
for L, S, J in ((1, 0.5, 1.5), (2, 1.0, 1.0)):
    for M in ([1.5, 0.5] if J == 1.5 else [1.0, 0.0]):
        s = sum(_l.clebsch_gordan(L, mL, S, mS, J, M)**2 * p2LM(L, mL)
                for mL in range(-L, L+1) for mS in np.arange(-S, S+0.5, 1.0))
        print("L=%d M=%+4.1f  <P2>=%+.6f  -5<P2>/Q_NN=%+.6f"
              % (L, M, s, -5*s/((3*M*M - J*(J+1))/3.0)))
print(b1_li7(), q_at(v7))
```

The full script as run, with all the tables of §5 formatted, is
`scratch/li7_rank2_estimate.py` in the working session; it is not committed,
because §9's excerpt regenerates everything the document quotes.

**The two numbers a reader should check first**, because everything else rests
on them:

```
L=1 M=+1.5  <P2>=-0.200000  -5<P2>/Q_NN=+1.000000      <- the 7Li coefficient
L=2 M=+1.0  <P2>=-0.100000  -5<P2>/Q_NN=+1.500000      <- the code's own 6Li DD coefficient,
                                                          b1_nuclear.cpp:356 writes 3*2/4 = 1.5
```

---

## 10. What this note does **not** claim

* It does **not** claim a b₁(⁷Li). It claims a leading estimate whose sign is
  not determined by the inputs in this tree (§5.4), and a second in-tree
  decomposition that disagrees with it (§5.5).
* It does **not** validate the light-cone convolution at A = 7. The quadrupole
  gate validates the **wave function**; nothing offline validates the
  convolution, and §4.1 explains why nothing can.
* ~~Q(⁷Li) = −4.00(3) fm² is **external and not in this tree**. It is quoted
  from memory of the standard compilations and **must be sourced before it is
  committed** (D11); the −4.06 fm² alternative is recorded beside it.~~
  **SOURCED AND COMMITTED 2026-09-06 (D11 paid, run 2026-09-06 task B3).** The
  constant is `LI7_QUADRUPOLE_FM2` = **−4.06 fm²** in `rc.hpp`, beside
  `LI6_QUADRUPOLE_FM2` — the −4.06 alternative, not the −4.00 this note quoted,
  because it comes from the **same evaluation** the ⁶Li constant does (TUNL's
  A = 5, 6, 7, Tilley *et al.*, NPA 708 (2002) 3, whose A = 7 half prints
  "Q = −40.6 ± 0.8 mb"), which is what "match conventions" means here. §4.2's
  ratio is therefore **0.858389** as the primary number and **0.871265**
  against the −4.00; the swap is 1.5 % and changes nothing this note concludes.
  See `../run_2026-09-06/phase_B_numbers.md` §B3 for the fetch, the eight-entry
  Stone compilation spread, and the measurement.
* The ⁶Li + n sketch of §5.5 is a **Clebsch–Gordan transfer only** — no orbital
  smearing, no struck neutron, from an r-space-only table. It is an
  order-of-magnitude cross-check, not a computation.
* No shipped number moved. The three gates reproduce the Phase D baseline
  exactly (top of this file).
