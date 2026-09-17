# Phase D2 — the structure-function injection map, and the design of
# `--unpol-sf` / `--pol-sf`

**Research only. No code, no test, no reference JSON, no CLI flag was changed by
this task.** Everything below is either read off the tree at commit `903fcc9`
(clean) or measured against the build in `build/`, with the recipe printed next
to the number.

Suites at the time of writing, all reproduced before this file was written:

| gate | result |
|---|---|
| `build/lipolgen_tests` | **390 cases / 17 217 416 assertions / 1 skipped / 0 failed** |
| `python -m pytest python/tests -q` | **345 passed** (49.86 s) |
| `python3 validation/check_physics_channels_links.py` | **1094 references checked (strict), 97 ranges, 6 external, 0 broken, 6 allow-listed** |

---

## 0. Verdict on the finding, and one correction to its premise

The physics half of D2 is **confirmed, in full**:

* Every structure function the shipped generator evaluates — F₂, F₁, F_L, g₁,
  g₂, R, and the EMC hook — is the toy backend on **every channel**:
  `ToyF2`, `ToyG1` (g₂ = Wandzura–Wilczek on that same g₁), `r_sigma_lt`, and
  an **empty** EMC hook.
* Phase A's `--b1-unpol {toy,ct18nlo,mstw}` reaches **exactly one** object,
  `Li6ConvolutionOptions::unpol` (`src/core/pipeline.cpp:1232` (anchor read 2026-09-15 "unpol = b1_unpol")), and nothing
  else. Verified numerically: with `--b1-model li6-convolution --b1-unpol
  {toy,mstw,ct18nlo}` the kernel's whole-nucleus F₂A is *bit-identical* to
  `NuclearF2(LI6(), ToyF2())` on every setting.
* The tagged channels have **no injection point at all**. `StruckClusterOptions`
  (`include/lipolgen/pipeline.hpp:207`) carries five fields —
  `inclusive_b1`, `delta_func`, `scenario`, `grid`, `with_perp` — and not one of
  them is an `UnpolSF`, a `PolSF` or an `RFunc`. `struck_cluster_kernel`
  (`src/core/pipeline.cpp:215` (as of adec442)) default-constructs `InclusiveKernel::Options`
  and sets only `b1_func` and `delta_func`. Neither the function nor
  `make_struck_cluster_source` is bound to Python at all.

**The premise that "the documents missed it" is false, and should not be
repeated.** `docs/PHYSICS_CHANNELS.md` §3 already states it, verbatim and
correctly, including the tagged clause:

> on the **tagged** channels `struck_cluster_kernel` exposes only
> `inclusive_b1` and `delta_func`, so F₂, g₁, R and EMC there are hard-locked to
> `ToyF2`/`ToyG1`/`r_sigma_lt`/none.

and the EMC row (`docs/PHYSICS_CHANNELS.md:193`) already says
`Options::emc_ratio` "is empty in every default kernel". `docs/USAGE.md`
§§2a and the `B1UnpolSource` block in `include/lipolgen/pipeline.hpp:395-406`
say the same about the scope of `--b1-unpol`. What is missing is not the
*documentation* of the gap — it is the *injection point*. D2 is a code task,
not a documentation task, and the phase-F rewrite of `PHYSICS_CHANNELS.md`
should carry the existing sentences forward, not "correct" them.

---

## 1. The map

### 1a. Every site that constructs or defaults an `InclusiveKernel`, `UnpolSF`, `PolSF`, `TensorSF` or R-function

`—` means the site does not set that slot at all, so the value is whatever the
row below it in the chain supplies. `(inherit)` means it takes the object from
another site rather than making one.

| # | site | builds | F₂ (`UnpolSF`) | g₁ (`PolSF`) | R (`RFunc`) | EMC hook | rank-2 (`TensorSF`) |
|---|---|---|---|---|---|---|---|
| **K1** | `src/core/xsec.cpp:57-85` `InclusiveKernel::InclusiveKernel(Ion, Options)` | the class-level defaults | `options.f2_source`, **null → a fresh `ToyF2`** (`:61-63`) | `options.g1_model`, **null → `ToyG1(nf2_.base(), options.r_func)`** (`:81`) | `options.r_func`, **null → `r_sigma_lt`** via `resolve_r` | `options.emc_ratio`, null → 1 | `b1_func`/`b2_func`/`delta_func`/`b*_32_func`/`b3`/`b4`, all null → 0 |
| **K2** | `src/core/pipeline.cpp:672-674` (as of adec442) `default_inclusive_kernel(const Ion&)` | forwards to K3 with `(Miller, 1.0, 1.0, nullptr)` | (inherit K3) | — | — | — | (inherit K3) |
| **K3** | `src/core/pipeline.cpp:676-738` `default_inclusive_kernel(ion, model, band, w_alpha_d, b1_unpol)` | **the inclusive and coherent kernel of every `Pipeline` run** | **hard-wired `std::make_shared<const ToyF2>()`** (`:762-763`) | **not set** → K1's `ToyG1` on that same `ToyF2` | **not set** → `r_sigma_lt` | **not set** → 1 | spin-1 only: `Li6B1(MillerB1)` / `Li6B1(CdksB1)`·band / `Li6ConvolutionB1`·band, plus `toy_delta_gluon(…,1e-2)` |
| **K4** | `src/core/pipeline.cpp:712-726` (as of adec442) the `Li6Convolution` branch inside K3 | `Li6ConvolutionOptions` | `o.unpol = b1_unpol ? b1_unpol : f2` (`:725` (as of adec442)) — **the one selectable slot in the whole tree** | n/a | **not set** → `Li6ConvolutionOptions::r_func` null → `r_sigma_lt` (`include/lipolgen/b1_nuclear.hpp:544`) | n/a | `deuteron_b1` not set → raw digitized CDKS column |
| **K5** | `src/core/pipeline.cpp:215-241` (as of adec442) `struck_cluster_kernel(channel, opt)` | **the DIS kernel of every tagged run** | **not set** → K1's `ToyF2` | **not set** → K1's `ToyG1` | **not set** → `r_sigma_lt` | **not set** → 1 | `b1_func = toy_b1` only when `opt.inclusive_b1`; `delta_func = opt.delta_func` |
| **K6** | `src/core/pipeline.cpp:951-956` the `Pipeline`'s non-tagged branch | picks `cfg_.kernel` if non-null, else K3 | (inherit) | (inherit) | (inherit) | (inherit) | (inherit) |
| **K7** | `src/core/pipeline.cpp:1541-1542` the `Pipeline`'s tagged branch | `make_struck_cluster_source` → K5. **`cfg_.kernel` is not read on this branch at all** | (inherit K5) | (inherit K5) | (inherit K5) | (inherit K5) | (inherit K5) |
| **K8** | `src/core/pipeline.cpp:935` (as of adec442) T1 breakup wiring | `BreakupOptions::f2` | **(inherit)** `dis_sampler_->kernel().nuclear_f2().base()` — correct by construction | n/a | n/a | n/a | n/a |
| **K9** | `src/core/breakup.cpp:127-129` `ClusterBreakup` ctor | its own fallback | `opt_.f2`, **null → a fresh `ToyF2`** | n/a | n/a | n/a | n/a |
| **K10** | `src/pythia/pythia_bridge.cpp:820-823` (as of 2e404b8) `PythiaBridge` ctor | the T2 struck-nucleon species draw | `opt.f2_source`, **null → a fresh `ToyF2`**; **nothing in the tree ever sets it** | n/a | n/a | n/a | n/a |
| **K11** | `src/core/b1_nuclear.cpp:479-484` `DeuteronConvolutionB1` ctor | the **A = 2 validation gate** object, not a pipeline kernel | `opt_.unpol`, null → `ToyF2` | n/a | `opt_.r_func`, **null → `r1998`** — deliberately different from K4 | n/a | n/a |
| **K12** | `src/core/b1_nuclear.cpp:591-600` `Li6ConvolutionB1` ctor | its own F₁ᵈ and F₁^α | `opt_.unpol`, null → `ToyF2` | n/a | `opt_.r_func` (null → `r_sigma_lt`) | n/a | `cdks_b1_raw_per_nucleon()` |
| **K13** | `examples/generate_inclusive.cpp:100-111` | the example's own kernel | not set → `ToyF2` | not set → `ToyG1` | not set | not set | `Li6B1(MillerB1)` + toy Δ for spin 1 |
| **K14** | `src/core/coherent.cpp` — **no site** | `CoherentSampler` (`src/core/coherent.cpp:322-328`) takes no structure function of any kind; `coherent.hpp` does not even include `sf.hpp` | n/a | n/a | n/a | n/a | n/a |

Sites that **consume** a backend without owning one (they inherit and therefore
need no new knob):

| site | what it reads |
|---|---|
| `src/core/sampler.cpp:148-166` `InclusiveSampler` cell cross sections | `kernel_->nuclear_f2()` and **its** `r_func()` |
| `src/core/generator.cpp:47-55` `InclusiveGenerator::proton_fraction` | `sampler_->kernel().nuclear_f2().base()` |
| `src/core/breakup.cpp:176-184` `ClusterBreakup::proton_fraction` | its own `f2_`, fed from K8 |
| `src/core/pipeline.cpp:979` (as of adec442) the coherent channel's cell weights | `dis_sampler_->cell_xsec_pb()` — i.e. **K3's F₂ and R**, reweighted by `f_coh(x)` |
| `src/core/rc.cpp:1140` (as of d3ac125) | `dis_->kernel().dsigma_unpol(x, q2, s)` |

### 1b. Reachability of each override

| slot | override exists in C++? | reachable from `PipelineConfig`? | from `lipolgen.make_config`? | from the CLI? |
|---|---|---|---|---|
| K3 `f2_source` (inclusive + coherent F₂/F₁/F_L) | yes, `InclusiveKernel::Options::f2_source` (`include/lipolgen/xsec.hpp:197`) | **only** by hand-building a whole kernel and assigning `PipelineConfig::kernel` (`include/lipolgen/pipeline.hpp:430` (as of adec442)) | **no** — `make_config` has no `kernel=` parameter | **no** |
| K3 `g1_model` | yes (`include/lipolgen/xsec.hpp:199`) | same — whole-kernel only | **no** | **no** |
| K3 `r_func` | yes (`include/lipolgen/xsec.hpp:231`) | same — whole-kernel only | **no** | **no** |
| K3 `emc_ratio` | yes (`include/lipolgen/xsec.hpp:230`) | same — whole-kernel only | **no** | **no** |
| K4 `Li6ConvolutionOptions::unpol` | **yes, by name** — `PipelineConfig::b1_unpol` + `b1_unpol_sf` (`include/lipolgen/pipeline.hpp:637,646`) | **yes** | **yes**, `b1_unpol=` (`python/lipolgen/__init__.py:436-447`) | **yes**, `--b1-unpol` (`python/lipolgen/cli.py:277` (as of adec442)) |
| K4 `Li6ConvolutionOptions::r_func` | yes, on the struct | **no** — `default_inclusive_kernel` never sets it | no | no |
| K5 tagged `f2_source` / `g1_model` / `r_func` / `emc_ratio` | **NO OVERRIDE EXISTS ANYWHERE.** `StruckClusterOptions` has no such field, and `struck_cluster_kernel` is not bound to Python | **no** | **no** | **no** |
| K5 tagged `b1_func` | only the boolean `inclusive_b1` (`toy_b1` or nothing) | yes | yes, `inclusive_b1=` | yes, `--inclusive-b1` |
| K5 tagged `delta_func` | yes, `StruckClusterOptions::delta_func` | yes, `cfg.struck` | no | no |
| K10 `PythiaBridgeOptions::f2_source` | yes (`include/lipolgen/pythia_bridge.hpp:284`) | n/a (the bridge is not in `PipelineConfig`) | **no — the field is not even bound**; `python/bindings.cpp:4076-4100` binds every other `PythiaBridgeOptions` field and omits this one | **no** |
| K11 / K12 `unpol` / `r_func` | yes, on the option structs, bound at `python/bindings.cpp:1497,1599` | n/a — gate objects, not pipeline kernels | n/a | n/a |

**The one-line summary of the map:** every row of table 1b is "no" except one.
The single structure-function backend the run surface can select is
`--b1-unpol` → K4's `Li6ConvolutionOptions::unpol`, and it is the one that
reaches the *smallest* observable — the b₁ numerator of the tensor weight,
deliberately not the rate. F₂, g₁, R and the EMC hook are unreachable on every
channel except through a hand-built `PipelineConfig::kernel`, which is
C++/pybind-only, is not offered by `make_config` or the CLI, and — see §7(i) —
is silently ignored on the tagged channels anyway.

---

## 2. What already inherits, and must not get a second knob

`docs/CONVENTIONS.md`'s "no physics number is defined twice" rule already does a
lot of the work here, and a design that adds a knob per site would break it. The
following are already correct and must be **left alone**:

* **K8** — the T1 breakup's species draw takes the sampler kernel's own
  `UnpolSF` object. It follows any new selector for free.
* **`InclusiveGenerator::proton_fraction`** — same, through the same accessor.
* **`InclusiveSampler`'s cell cross sections and R** — read off the kernel's
  `NuclearF2`, so both move together with the kernel.
* **The coherent channel** — has no structure function of its own; it rides
  K3's cell cross sections. A selector on K3 reaches it automatically, and
  `CoherentSampler` must not grow an `UnpolSF` argument.
* **g₂** — `PolSF::g2p/g2n` is Wandzura–Wilczek built on *this backend's own*
  g₁ (`include/lipolgen/sf.hpp:175-177`), so a `--pol-sf` selector moves g₁ and
  g₂ consistently with one field. No `--g2-sf` is needed or wanted.

The one site that does **not** inherit and would become inconsistent the moment
a selector lands is **K10**, the PYTHIA bridge (see §5.6).

---

## 3. What the toy defaults cost, measured

All rows below are measurements made against `build/` on 2026-09-04 with
`source env.sh`. Every recipe is a few lines of Python through the installed
bindings; none needs a file that is not in the tree.

### 3a. The unpolarised rate

Sum of `InclusiveSampler::cell_xsec_pb()` over the 3051 accepted cells of the
shipped ⁶Li inclusive window (`generator_scenario(Scenario())`, default
`SamplerGrid`, `default_configs("6Li")[1]`), with `InclusiveKernel::Options::f2_source`
swapped and nothing else changed:

| `f2_source` | Σ cell σ [pb] | ratio to toy |
|---|---|---|
| `ToyF2` (**shipped**) | 591 846.2 | 1.0000 |
| `LhapdfSF("CT18NLO", 0)` | 472 571.9 | **0.7985** |
| `MstwSF()` | 469 556.5 | **0.7934** |

The same swap on the **tagged** struck-cluster kernel (deuteron target, the ⁶Li
beam's per-nucleon momentum) gives *the identical* 0.7985, because both ⁶Li and
the embedded deuteron are isoscalar and every structure function here is
per nucleon. The triton and free-nucleon `dis_target`s are not isoscalar and
were **not** measured.

Recipe:

```python
from lipolgen import _lipolgen as _l
beams = _l.default_configs('6Li')[1]; sc = _l.generator_scenario(_l.Scenario())
for f2 in (None, _l.LhapdfSF('CT18NLO', 0), _l.MstwSF()):
    o = _l.InclusiveKernel.Options()
    if f2 is not None: o.f2_source = f2
    s = _l.InclusiveSampler(_l.InclusiveKernel(_l.li6(), o), beams, sc, _l.SamplerGrid())
    print(sum(s.cell_xsec_pb))
```

Pointwise F₂ᵖ at Q² = 10, toy → CT18NLO / MSTW2008 LO:

| x | ToyF2 | CT18NLO (ratio) | MSTW08LO (ratio) |
|---|---|---|---|
| 0.01 | 0.66627 | 0.61397 (0.9215) | 0.5433 (0.8154) |
| 0.05 | 0.4673 | 0.47244 (1.0110) | 0.42386 (0.9070) |
| 0.1 | 0.38524 | 0.42742 (1.1095) | 0.39035 (1.0133) |
| 0.2 | 0.2844 | 0.36621 (1.2877) | 0.34462 (1.2118) |
| 0.3 | 0.20795 | 0.28508 (**1.3709**) | 0.28382 (1.3649) |
| 0.5 | 0.091503 | 0.11566 (1.2641) | 0.13249 (**1.4479**) |
| 0.7 | 0.023348 | 0.022509 (0.9641) | 0.030179 (1.2926) |

and F₂ⁿ/F₂ᵖ: toy 0.9625 / 0.8500 / 0.6250 at x = 0.05 / 0.2 / 0.5 against
CT18NLO 0.9219 / 0.7220 / 0.5035 and MSTW 0.9271 / 0.7427 / 0.4856. The toy
n/p ratio is a straight line `clip(1 − 0.75x, 0.25, 1)`; it is up to **24 %**
above the fits at x = 0.5, and it is what the T1 and T2 species draws use.

### 3b. The polarised sector

`InclusiveKernel::Options::g1_model` swapped from the default `ToyG1` to
`LhapdfG1("NNPDFpol11_100", 0)`, ⁶Li, per-nucleon g₁A and A_∥ at y = 0.5 with
the default `target_mass = true`:

| x | Q² | g₁A toy | g₁A NNPDFpol | ratio | A_∥ ratio |
|---|---|---|---|---|---|
| 0.01 | 2.5 | −0.0847927 | −0.0567341 | 0.6691 | 0.6691 |
| 0.05 | 5 | 0.0323251 | 0.0218549 | 0.6761 | 0.6756 |
| 0.1 | 10 | 0.0343653 | 0.0382812 | 1.1139 | 1.1137 |
| 0.2 | 10 | 0.0255593 | 0.0344831 | **1.3491** | 1.3493 |
| 0.3 | 15 | 0.0184751 | 0.0239836 | 1.2982 | 1.2986 |
| 0.5 | 25 | 0.00788548 | 0.00773991 | 0.9815 | 0.9818 |
| 0.7 | 50 | 0.00196416 | 0.00124874 | **0.6358** | 0.6358 |

So the polarised swap is worth **×0.64 to ×1.35** on the ⁶Li longitudinal
asymmetry over the generator window — not a normalisation, and non-monotone.

The **neutron** is worse than the nucleus, and the reason is a sign. At
Q² = 10:

| x | ToyG1 g₁ⁿ | NNPDFpol g₁ⁿ | ratio |
|---|---|---|---|
| 0.05 | −0.242789 | −0.248012 | 1.0215 |
| 0.1 | −0.0800273 | −0.136452 | 1.7051 |
| 0.2 | −0.0113549 | −0.0607773 | **5.3525** |
| 0.3 | **+0.00520678** | **−0.0273959** | **−5.2616** |
| 0.5 | +0.00778817 | −0.000370697 | −0.0476 |
| 0.7 | +0.00247092 | +0.00221761 | 0.8975 |

`ToyG1`'s `a1n(x) = −0.07(1−x)² + 0.8x^2.2` crosses zero and goes positive
around x ≈ 0.25; NNPDFpol1.1's g₁ⁿ stays negative to x ≈ 0.6. **The shipped
toy g₁ⁿ has the wrong sign over roughly 0.25 < x < 0.6.** That matters for
every neutron-tagged observable (`d`+`p` tagging, whose `dis_target` is
`NEUTRON_TARGET`) far more than it does for the isoscalar ⁶Li, where the
proton term dominates the sum. This is a **first measurement in this run** and
it is the single strongest argument for D2 being worth doing at all.

### 3c. `--unpol-sf` alone already moves g₁ — and only partly

Because `InclusiveKernel` builds its default `ToyG1` **on its own base
`UnpolSF`** (`src/core/xsec.cpp:81`), swapping `f2_source` alone moves F₁, g₁
and A_∥ together. With `f2_source = CT18NLO` and `g1_model` left at the
default, ⁶Li at y = 0.5:

| x | Q² | F₁ ratio | g₁ ratio | A₁ = g₁/F₁ ratio | A_∥ ratio |
|---|---|---|---|---|---|
| 0.05 | 5 | 0.9942 | 1.0544 | 1.0605 | 1.0603 |
| 0.1 | 10 | 1.0711 | 1.1346 | 1.0593 | 1.0591 |
| 0.3 | 15 | 1.2038 | 1.3050 | 1.0841 | 1.0843 |
| 0.5 | 25 | 1.0014 | 1.0551 | 1.0537 | 1.0540 |

A₁ does **not** stay fixed (it moves 5–8 %) because `ToyG1` combines
`a1p(x)·F₁ᵖ` and `a1n(x)·F₁ⁿ` with the backend's own n/p split, which differs
between fits. This is a design constraint, not a bug: **the two selectors are
not orthogonal**, and any documentation of them must say which one a plot used
*and* what the other was set to.

### 3d. R

Nothing in the pipeline ever sets an `RFunc`, so R is `r_sigma_lt` at all four
places `InclusiveKernel` needs it. Over the 3051 accepted cells of the shipped
⁶Li run:

* `r_sigma_lt` spans **0.0046 – 0.1763**; `r1998` spans **0.0124 – 0.4010**.
* `(1 + r1998)/(1 + r_sigma_lt)` spans **0.9203 – 1.1910**, cell-σ-weighted mean
  **1.1229**.
* **38.18 %** of the accepted cell cross section lies **outside** R1998's own
  stated support (`R1998_X_MIN/MAX` = 0.005/0.86, `R1998_Q2_MIN/MAX` = 0.5/130,
  `include/lipolgen/sf.hpp:52-55`), where `r1998(..., clip = true)` returns the
  clipped boundary value.

So R1998 is **not** a drop-in replacement default: over more than a third of
the window it is an extrapolation of a fit past its own stated range. R must
stay its own axis (§5.7).

### 3e. The LHAPDF grid boundary — a hard constraint on `--unpol-sf ct18nlo`

| set | XMin | Q range | fraction of the shipped ⁶Li cell σ below the grid |
|---|---|---|---|
| CT18NLO | 1e−9 | Q ≥ 1.295 GeV, i.e. **Q² ≥ 1.677** | **42.33 %** |
| NNPDFpol11_100 | 1e−5 | Q ≥ 1 GeV (Q² ≥ 1), Q ≤ 316 GeV | **0.00 %** |
| EPPS21nlo_CT18Anlo_Li6 | 1e−7 | Q ≥ 1.3 GeV | (not used by a kernel) |

The accepted cells run Q² = 1.054 … 1897 and x = 2.884e−4 … 0.955. The
scenario cut is `q2_min = 0.7`; the *grid* is what stops at 1.054.

Below its grid CT18NLO does **not** freeze — LHAPDF continues the evolution
downward, and it falls fast. F₂ᵖ at x = 3e−4: toy 0.6812 / CT18NLO 0.2916 at
Q² = 0.7 (a factor **2.34** low), 0.6812 / 0.4008 at Q² = 1.0. MSTW is milder
(0.4923 and 0.5499 at the same two points) but is also below its grid.

**Consequence for the design:** `--unpol-sf ct18nlo` puts 42 % of the shipped
window's rate on a silent below-grid extrapolation. That must be *said* — at
the run banner and in the npz meta — or the flag must refuse to run without an
explicit `--q2-min` raise. It must never be silent, which is the same rule
`--b1-unpol` already enforces for a missing tier.

---

## 4. How Phase A threaded `B1UnpolSource` — the pattern, read off the tree

Seven pieces, and the generalisation should copy all seven:

1. **A `std::uint8_t` enum with a `Custom` value**, declared next to the field
   it governs, with the whole "why this exists / what each value is / what the
   scope is and what it costs" block in the header comment
   (`include/lipolgen/pipeline.hpp:350-408`). `Custom` is the escape hatch and
   the provenance label at once, copied from `OpticsChoice`.
2. **Two fields, kept in step**: the enum (`b1_unpol`, the PROVENANCE — it is
   what `meta` records) and the object (`b1_unpol_sf`, the REALISATION)
   (`include/lipolgen/pipeline.hpp:637,646`).
3. **The core library cannot build the named backends and does not try.**
   `MstwSF` is in the PYTHIA tier, `LhapdfSF` in the LHAPDF tier, and `sf.hpp`'s
   rule is that the core links neither. So the object is built one layer up, in
   `_lipolgen.set_b1_unpol(cfg, source)` (`python/bindings.cpp:4150-4188`),
   which is the `set_pythia_hadronizer` arrangement. A missing tier throws
   there with the tier named.
4. **`validate()` refuses two things** (`src/core/pipeline.cpp:598-668`):
   a named backend the flag would not reach ("a knob that did not run may not be
   recorded as if it had"), and a named backend with an empty object slot
   ("never a silent fallback to `ToyF2`"). Both messages name the flag, the
   reason, and the fix.
5. **`meta` is written unconditionally**, with a `caller-supplied kernel`
   string when `cfg.kernel` wins (`python/bindings.cpp:479-481` (as of adec442)). The
   unconditional write is *why* rule 4's first clause exists.
6. **The Python surface keeps the pair consistent in both directions**:
   assigning an object sets `Custom`, clearing it sets `Toy`
   (`python/bindings.cpp:4580-4595`); `set_b1_unpol` goes the other way.
7. **The CLI checks the tier at the command line** before the binding can throw
   three frames down (`require_b1_unpol_tier`, `python/lipolgen/cli.py:448-475`),
   and the **run banner prints the setting on every run, the default included**
   (`python/lipolgen/cli.py:1588-1599`) — "a run that does not say which one it
   used is not reproducible from its own log".

---

## 5. The design

### 5.1 The two new selectors

```cpp
// include/lipolgen/pipeline.hpp, next to B1UnpolSource
enum class UnpolSfSource : std::uint8_t { Toy, Mstw, Ct18Nlo, Custom };
const char* unpol_sf_name(UnpolSfSource s);

enum class PolSfSource   : std::uint8_t { Toy, NnpdfPol, Custom };
const char* pol_sf_name(PolSfSource s);

// ... in PipelineConfig
UnpolSfSource unpol_sf = UnpolSfSource::Toy;
std::shared_ptr<const UnpolSF> unpol_sf_obj;     // empty for Toy
PolSfSource   pol_sf   = PolSfSource::Toy;
std::shared_ptr<const PolSF>  pol_sf_obj;        // empty for Toy
```

`Toy` on both leaves the corresponding `Options` slot **null**, which is
exactly what K3 and K5 do today, so the default path is byte-identical by
construction rather than by inspection.

**Polarised sets available.** `ls ../deps/install/share/LHAPDF/` gives
`CT18NLO`, `EPPS21nlo_CT18Anlo_Li6`, `NNPDFpol11_100`, plus `lhapdf.conf` and
`pdfsets.index`. **`NNPDFpol11_100` is the only polarised set installed**, and
it is already the declared default of `LhapdfG1`
(`include/lipolgen/sf.hpp:218` (as of ac22331)). So `--pol-sf {toy,nnpdfpol}` is the complete
list on this machine and there is no second row to offer. **Note that
`CT18ANLO` is not installed** — see §6.

### 5.2 Which kernels must change

| kernel | change? | why |
|---|---|---|
| **K3** `default_inclusive_kernel(ion, model, band, w, b1_unpol)` | **YES** — two new trailing arguments, `unpol_sf` and `pol_sf`, both defaulting to `nullptr` | it is the inclusive **and** coherent kernel; the coherent channel inherits for free through `dis_sampler_->cell_xsec_pb()` |
| **K5** `struck_cluster_kernel(channel, opt)` | **YES** — three new `StruckClusterOptions` fields: `f2_source`, `g1_model`, and (see §5.7) `r_func`; `Pipeline` fills them from the config at `src/core/pipeline.cpp:850-853` (as of adec442) beside `sopt.scenario` and `sopt.grid` | this is the gap D2 names; it is the *only* way a tagged run can reach a real backend |
| **K2** `default_inclusive_kernel(const Ion&)` | **NO** — leave the one-argument overload exactly as it is | it is pinned bit-for-bit by `validation/reference/b1_default_li6.json` at rtol 1e-12 through `tests/test_b1_nuclear.cpp:1851` (T9) and `python/tests/test_b1_model.py:339`. Adding defaulted arguments to K3 keeps this call site untouched, which is how Phase A already did it |
| **K1** `InclusiveKernel`'s own defaults (`src/core/xsec.cpp:57-85`) | **NO — absolutely not** | `validation/reference/xsec.json` is dumped from a *default-constructed* Python kernel and rebuilt in C++ at `tests/test_reference.cpp:48-82` with `f2_source`, `g1_model` and `r_func` all left at their defaults. Any change to the null-branch of `:62-64` (anchor read 2026-09-15) or `:97` (anchor read 2026-09-15 "return -1e-2") moves that gate. The selectors must live one layer up and hand the class an explicit object; they must never change what "null" means |
| **K9** `ClusterBreakup`'s `ToyF2` fallback | **NO** | `Pipeline` already overwrites it from the kernel (`src/core/pipeline.cpp:1687`). The fallback only fires for a direct caller, and a second knob here would be a second definition of the same number |
| **K11** `DeuteronConvolutionB1` (the A = 2 gate) | **NO** | its `unpol` default is `ToyF2` and its `r_func` default is `r1998` **on purpose** — it exists to reproduce CDKS's own Fig. 4, and its inputs are CDKS's choices. Sweeping it into a general selector would silently move the G3a/G3b verdict rows that `docs/USAGE.md` §2a and `phase_A_numbers.md` quote |
| **K13** `examples/generate_inclusive.cpp` | **NO** | it is the documented "same choice as `default_inclusive_kernel`" example and is not a run surface |
| **K14** `CoherentSampler` | **NO** | it has no structure function; it must not grow one |

### 5.3 What stays bit-for-bit, and how it is proved rather than asserted

* **The default run.** `unpol_sf = Toy` and `pol_sf = Toy` leave both `Options`
  slots null on K3 and K5, so `InclusiveKernel` takes exactly the branches it
  takes today. Nothing in `src/core/xsec.cpp` changes.
* **The rtol 1e-12 reference gates.** All eight files in
  `validation/reference/_manifest.json` — `xsec.json`, `spin.json`,
  `beams.json`, `tagged.json`, `spectator.json`, `coherent.json`,
  `bookkeeping.json`, `b1_default_li6.json` — are compared against kernels
  built either from K1's raw defaults (`tests/test_reference.cpp:48-82`, which
  never touches `f2_source`, `g1_model` or `r_func`) or from K2's one-argument
  overload (`tests/test_b1_nuclear.cpp:1864` (anchor read 2026-09-15 "const auto kernel")). **No test in the tree calls
  `struck_cluster_kernel` directly** (`grep` finds it only in three comments at
  `tests/test_b1_nuclear.cpp:1899-1905`), and none reaches K3's five-argument
  form except through the `b1_unpol` tests that already exist. So no reference
  file can move, `_manifest.json` needs no regeneration, and `--record-ranges`
  is not involved.
* **The tagged pipeline tests DO exercise K5** — `tests/test_tagged.cpp`,
  `tests/test_pipeline.cpp` and `python/tests/test_pipeline.py` all run tagged
  channels end to end. They are protected not by a reference JSON but by the
  null-default argument above, and the gate that turns that from an assertion
  into a proof is the first pytest of §5.10.
* **The proof obligation** is the same one Phase A carried: a pytest that runs
  the pipeline with the flag absent and with the flag explicitly at `toy`, and
  asserts `np.array_equal` on `cell_xsec_pb`, `sigma_per_category_pb()` and the
  generated `x`/`q2`/`y`/`weight` columns. That test already exists for
  `--b1-unpol` at `python/tests/test_b1_model.py:693-704`; copy it.

### 5.4 `b1_unpol`: keep it separate, and add one refusal

**Recommendation: keep `--b1-unpol` as a separate flag. Do not alias it.**

The reasons are physics reasons, not compatibility ones:

1. **They are different quantities with different conventions.**
   `--unpol-sf` sets the whole-nucleus F₂A/F₁A of the *rate*, through
   `NuclearF2` and the massless `f1a`. `--b1-unpol` sets the *deuteron's*
   per-nucleon F₁ inside CDKS Eq. (22), through `f1_cdks`, which is explicitly
   **not** `NuclearF2::f1a` (`include/lipolgen/b1_nuclear.hpp:500-503`,
   design 2.4). Two objects, two conventions, at two different Q² regimes.
2. **`--b1-unpol` is a CDKS-comparability choice, and it is load-bearing.**
   The A = 2 gate passes at G3b = 0.843243 with MSTW2008 LO and fails at 0.440
   on the toy. That row is a statement about reproducing a published figure. A
   user who wants a realistic *rate* (`--unpol-sf ct18nlo`) must not be forced
   to move the b₁ gate row off MSTW as a side effect, and vice versa.
3. **The tree already says so, in nine places.**
   `include/lipolgen/pipeline.hpp:426-437` documents the scope of `--b1-unpol`
   as "b₁ and NOTHING else" and lists the price. Aliasing would invalidate that
   paragraph and the **eight** other sites that repeat the same sentence. The
   full list, assembled by grepping `f2_source`, `b1_unpol` and
   `B1UnpolSource` across the tree (excluding `build/` and `.skbuild/`) and
   keeping the hits that assert the scope:

   | # | site | form |
   |---|---|---|
   | 1 | `include/lipolgen/pipeline.hpp:486-497` | the `B1UnpolSource` header block (the origin) |
   | 2 | `src/core/pipeline.cpp:1211-1216` (anchor read 2026-09-15) | the comment in the `Li6Convolution` branch |
   | 3 | `python/bindings.cpp:4079` | the `B1UnpolSource` enum docstring |
   | 4 | `python/lipolgen/__init__.py:203` (anchor read 2026-09-15 "B1_UNPOL = {") | the `B1_UNPOL` comment |
   | 5 | `python/lipolgen/cli.py:547` (anchor read 2026-09-16 `add_argument("--b1-unpol"`) | the `--b1-unpol` help string |
   | 6 | `python/lipolgen/cli.py:1401-1403` (anchor read 2026-09-15) | the run banner |
   | 7 | `python/tests/test_b1_model.py:838` (anchor read 2026-09-15 "assert tm.f1") | the assertion `tm.f1 == tt.f1` and its comment |
   | 8 | `tests/test_b1_nuclear.cpp:2082` (anchor read 2026-09-15 "ta.f1 == tb.f1") | the same assertion in C++ |
   | 9 | `docs/USAGE.md:562-568` | "Scope, and the price of it" |

   `docs/PHYSICS_CHANNELS.md:195-195` (anchor read 2026-09-15 "as `LhapdfSF") states the *consequence* rather than the
   sentence, so it is a tenth site to check by hand and no single grep pattern
   catches all ten — which is itself the reason to keep the two flags separate
   rather than re-derive this paragraph.

**But there is a collision that must be closed.** Today
`src/core/pipeline.cpp:1232` (anchor read 2026-09-15 "move(b1_unpol") reads

```cpp
o.unpol = b1_unpol ? std::move(b1_unpol) : f2;   // f2 is the kernel's own ToyF2
```

If `--unpol-sf` replaces `f2` with `CT18NLO`, then `b1_unpol = Toy` silently
starts meaning "CT18NLO", while `meta["b1_unpol"]` still writes `"toy"`
(`python/bindings.cpp:479-481` (as of adec442)). **That is precisely the defect the whole
provenance discipline exists to prevent** — a run labelled `toy` whose b₁ is
the CT18NLO one.

Two ways out, in order of preference:

* **(a) Refuse the ambiguity.** `validate()` throws when
  `b1_model == Li6Convolution && unpol_sf != UnpolSfSource::Toy &&
  b1_unpol == B1UnpolSource::Toy`, with a message that says: the kernel's
  `f2_source` is no longer `ToyF2`, so `b1_unpol = toy` no longer names what
  the convolution will fold against — say it explicitly. This costs one
  refusal, changes no existing meaning, keeps the default meta string `"toy"`
  stable, and is the same shape as the four refusals already in that function.
  **This is the recommendation.**
* **(b) Add a `B1UnpolSource::Inherit` value** meaning "the kernel's
  `f2_source` object, shared", and make it the new default. Truthful and
  arguably cleaner, but it changes `meta["b1_unpol"]` on a *default* run from
  `"toy"` to `"inherit"`, which breaks
  `python/tests/test_b1_model.py`'s meta assertions and makes every existing
  npz's meta string ambiguous against new ones. Not worth it.

### 5.5 What `meta` records

Three new keys, written **unconditionally** beside the existing b₁ block at
`python/bindings.cpp:463-482`, with the same `caller-supplied kernel` guard:

```cpp
meta["unpol_sf"] = p.config().kernel ? std::string("caller-supplied kernel")
                                     : std::string(unpol_sf_name(p.config().unpol_sf));
meta["pol_sf"]   = p.config().kernel ? std::string("caller-supplied kernel")
                                     : std::string(pol_sf_name(p.config().pol_sf));
```

and, because §3e makes it load-bearing, one boolean:

```cpp
// true when the run's accepted (x, Q2) window extends below the selected
// grid's own QMin^2 -- 42.33 % of the shipped 6Li window for CT18NLO.
meta["unpol_sf_below_grid"] = <computed from the sampler's q2_cells>;
```

The reason is the one already written at `:474-478`: two otherwise identical
npz files must not be indistinguishable in `meta`. Here the difference is a
factor 0.80 on the total cross section, so the key is not optional.

`meta["b1_unpol"]` keeps its current meaning and its current
unconditional write.

### 5.6 The PYTHIA bridge — the one site that breaks the moment this lands

`PythiaBridgeOptions::f2_source` (`include/lipolgen/pythia_bridge.hpp:284`)
selects the T2 struck-nucleon species draw, and:

* **nothing in the tree ever sets it** — it is null on every path, so the draw
  is on `ToyF2` (`src/pythia/pythia_bridge.cpp:863-866`);
* **it is not bound to Python** — `python/bindings.cpp:4995-5017` binds
  `seed`, `settings`, `verbosity`, `max_retries`, `headroom`, the three
  `include_*`, `q2_pdf_min`, `with_neutron_instance`, `coherent_t2`, `pom_set`,
  `pom_rescale` and `nucleon_choice`, and omits `f2_source` alone;
* the CLI builds the bridge itself at `python/lipolgen/cli.py:545-552` (as of a7b3d18), so
  `Pipeline` never sees it and cannot forward the kernel's backend the way it
  forwards it to the breakup at `src/core/pipeline.cpp:1687`.

Today that is harmless, because both ends are `ToyF2`. **The moment
`--unpol-sf ct18nlo` exists it is a live inconsistency**: the T0 rate and the
T1 species draw would be on CT18NLO while the T2 species draw stays on the toy,
and the toy's F₂ⁿ/F₂ᵖ is up to 24 % away from CT18NLO's at x = 0.5. The header
of `f2_source` already asks for the opposite ("hand it the same backend the
kernel was built with to keep the two draws consistent").

**So `--unpol-sf` must land with three extra changes:** bind
`PythiaBridgeOptions::f2_source`, set it in `cli.py` from the same selector
right where `popts.pom_set` is set, and add a pytest that a
`--hadronize --unpol-sf ct18nlo` run has the bridge and the kernel on the same
object. Without them the flag ships a known self-inconsistency.

### 5.7 R, and what this design deliberately does **not** cover

R is a **third** axis and is not part of `--unpol-sf` or `--pol-sf`. §3d
measures why it cannot simply follow them: `r1998` differs from `r_sigma_lt` by
up to 19 % in `(1+R)` and 38 % of the shipped window's rate is outside R1998's
own support. `Li6ConvolutionOptions::r_func`'s header already states the clean
fix and states that it was not made
(`include/lipolgen/b1_nuclear.hpp:505-506` (as of adec442)):

> The clean fix, NOT made here, is to thread ONE R hook through
> `default_inclusive_kernel` into both this and the kernel's `UnpolSF`.

Recommendation: add `RFunc` plumbing to K3 and K5 **as a nullable field with no
CLI flag yet** (so the C++ caller and a future `--r-model` have one place to
put it), and leave the selector itself to a separate task with its own measured
justification. Do not ship an `--r-model r1998` in the same change as
`--unpol-sf`: two independent 10–40 % movements landing together are
unattributable.

Also **not** covered, and deliberately:

* the **EMC hook** (`Options::emc_ratio`) — see §6;
* the polarised-EMC `medium_ratio` argument of `PolSF::g1_nucleus`, which
  `InclusiveKernel::g1a` never passes (`src/core/xsec.cpp:88`);
* **7Li's empty rank-2 slots** — that is D1, and `--pol-sf` does not touch it.

### 5.8 What `validate()` must refuse

Five clauses, all modelled on `src/core/pipeline.cpp:450-520` (anchor read 2026-09-15):

1. **Named-but-empty.** `unpol_sf != Toy && !unpol_sf_obj` → throw, naming the
   missing tier (PYTHIA for `mstw`, LHAPDF for `ct18nlo`) exactly as
   `:640-667` does. Same for `pol_sf != Toy && !pol_sf_obj` (LHAPDF for
   `nnpdfpol`). **Never a silent fallback to the toy.**
2. **Toy-with-an-object.** `unpol_sf == Toy && unpol_sf_obj` → throw; assigning
   an object from Python must set `Custom` (the `:3800-3807` setter pattern),
   so reaching this state means the two fields drifted.
3. **Caller-supplied kernel wins.** `kernel && (unpol_sf != Toy || pol_sf != Toy)`
   → throw. The kernel silently wins over the flags and `meta` would print both,
   which is the `:555-560` argument verbatim.
4. **The b₁ collision of §5.4(a).** `b1_model == Li6Convolution &&
   unpol_sf != Toy && b1_unpol == Toy` → throw, because `b1_unpol = toy` would
   name an object that is no longer `ToyF2`.
5. **The below-grid clause.** `unpol_sf == Ct18Nlo` with a scenario whose
   accepted window reaches below Q² = 1.677 is **not** refused — 42 % of the
   default window is there and refusing would make the flag unusable — but it
   **must** be recorded in `meta` (§5.5) and printed at the banner (§5.9). If
   the reviewer prefers a refusal, it belongs behind an explicit
   `--allow-below-grid`, never silent.

Note what must **not** be refused: `unpol_sf`/`pol_sf` on a **tagged** channel.
Unlike `b1_model`, these reach the rate on every channel, which is the whole
point of the task, so the "inclusive only" guard of `:561-570` has no analogue
here.

### 5.9 Banner and CLI

* `--unpol-sf {toy,ct18nlo,mstw}` and `--pol-sf {toy,nnpdfpol}` with
  `default=None` and `choices=sorted(...)`, beside `--b1-unpol` in
  `python/lipolgen/cli.py`.
* A `require_unpol_sf_tier` / `require_pol_sf_tier` pair copying
  `require_b1_unpol_tier` (`python/lipolgen/cli.py:448-475`) so the tier error
  arrives at the command line, not three frames down.
* Two `UNPOL_SF` / `POL_SF` dicts in `python/lipolgen/__init__.py` beside
  `B1_UNPOL`, and `unpol_sf=` / `pol_sf=` parameters on `make_config`
  (`python/lipolgen/__init__.py:429-430` and `:436-447`).
* `set_unpol_sf(cfg, source)` / `set_pol_sf(cfg, source)` bindings, copying
  `set_b1_unpol` (`python/bindings.cpp:4150-4188`) — the core cannot build
  either backend and must not try.
* **The banner prints both on every inclusive, coherent and tagged run,
  including the default**, with the measured factor and the below-grid warning.
  This is the `cli.py:576-587` rule: a run that does not say which backend it
  used is not reproducible from its own log. The banner must also say the
  §3c fact — that `--unpol-sf` alone moves g₁ through the default `ToyG1`, so
  `--pol-sf toy` is not "g₁ unchanged".

### 5.10 The pytest gates

Ten tests, all with a direct ancestor in `python/tests/test_b1_model.py`:

| # | test | model |
|---|---|---|
| 1 | absent flag == explicit `toy`: `cell_xsec_pb`, `sigma_per_category_pb()` and the `x`/`q2`/`y`/`weight` columns all `np.array_equal`, on **inclusive, coherent and one tagged channel** | `:693-704` |
| 2 | the two fields stay in step: assigning an object sets `Custom`, clearing it sets `Toy`, `set_*` builds the object, `Custom` cannot be asked for by name | `:707-728` |
| 3 | `validate()` refuses a named backend with an empty slot, for each of `mstw`/`ct18nlo`/`nnpdfpol`/`Custom`, with the name in the message | `:731-755` |
| 4 | `validate()` refuses `toy` with an object attached | `:750-755` |
| 5 | `validate()` refuses a caller-supplied kernel together with either selector | `:774-791` |
| 6 | `validate()` refuses the §5.4(a) b₁ collision | new |
| 7 | CLI: both flags in `DEFAULTS`, default `"toy"`, every choice named in `--help` | `:794-802` |
| 8 | `@needs_lhapdf`: the **tagged** struck-cluster kernel actually moves — the gate that D2 exists for. Assert `p.dis_sampler.kernel.nuclear_f2.f2a(0.3, 10)` differs from the `ToyF2` value and equals the `CT18NLO` one | new |
| 9 | `@needs_lhapdf`: `meta["unpol_sf"]`, `meta["pol_sf"]` and `meta["unpol_sf_below_grid"]` are written on every run, and read `"caller-supplied kernel"` when `cfg.kernel` is set | `:855-860` |
| 10 | `@needs_pythia @needs_lhapdf`: with `--hadronize --unpol-sf ct18nlo` the bridge's `f2_source` **is** the kernel's object (§5.6) | new |

Plus one C++ doctest asserting that `default_inclusive_kernel(LI6())` — the
one-argument overload — is still byte-identical to
`validation/reference/b1_default_li6.json`, which T9
(`tests/test_b1_nuclear.cpp:1851`) already does and which must be re-run rather
than re-recorded.

---

## 6. The EMC baseline — does a general selector interact with it?

**Mechanically, no. Conceptually, yes, and the trap is real.**

Mechanically:

* `EmcBaseline::Epps21` is the library default (`include/lipolgen/sf.hpp:309`)
  but it is realised by a **transcribed compile-time constant**,
  `EMC_VALENCE_DEPLETION_EPPS21 = 0.031052077003862335`
  (`include/lipolgen/sf.hpp:324`), returned by `emc_valence_depletion`
  (`src/core/sf.cpp:271-278`). It reads no `UnpolSF` and calls no LHAPDF grid at
  run time. A selector cannot move it.
* The functions that consume it — `cbt_valence_scale`, `tmt_valence_scale`,
  `cbt_polarized_emc_ratio`, `tmt_polarized_emc_ratio` — take no `UnpolSF`
  argument at all.
* **And none of them is wired to any kernel.** `InclusiveKernel::Options::emc_ratio`
  is set by *no* site in `src/`, `python/` or `examples/` — grep confirms the
  only occurrences are the declaration (`include/lipolgen/xsec.hpp:230`), the
  consumption (`src/core/xsec.cpp:65`) and the pybind property
  (`python/bindings.cpp:1535` (anchor read 2026-09-15 `arg("emc_ratio`)). So **no shipped run applies any EMC ratio to
  F₂A on any channel**, exactly as `docs/PHYSICS_CHANNELS.md:175` states.

Conceptually, the interaction is a **denominator-matching** trap that
`include/lipolgen/sf.hpp:309-313` already documents for a different reason:

> The denominator is EPPS21's OWN proton baseline, so the fit cancels and the
> ratio is the nuclear modification alone; it was CT18NLO until 2026-08-29,
> which mixed in the CT18A-vs-CT18 difference between two proton fits and made
> the depletion 4.2 % SHALLOWER (0.02979).

The EPPS21 depletion 0.031052… is `F2(EPPS21nlo_CT18Anlo_Li6) / F2(CT18ANLO)`.
Its denominator is **CT18ANLO**, not CT18NLO
(`include/lipolgen/sf.hpp:330`). If a future task wires an
`emc_ratio` hook — via `Epps21Ratio` (`include/lipolgen/lhapdf_sf.hpp:48-53`),
whose own `proton_set` default is **`"CT18NLO"`**, deliberately not CT18ANLO —
and the run is also on `--unpol-sf ct18nlo`, the ratio and the baseline would
be referenced to two different proton fits, which is the *same* 4.2 %-class
error the constant's provenance note describes.

Two concrete consequences for this design:

1. **`--unpol-sf` must not offer a `ct18anlo` row on this machine.**
   `ls ../deps/install/share/LHAPDF/` lists only `CT18NLO`,
   `EPPS21nlo_CT18Anlo_Li6` and `NNPDFpol11_100`. **CT18ANLO is not
   installed**, so the one value that would make an EPPS21 EMC hook
   self-consistent cannot be selected. Offering it would produce an
   `LHAPDF::readPDF` failure, not a fallback — which is the right failure mode,
   but it should not be advertised as a choice until the set is on disk.
2. **The `--unpol-sf` header block must carry a one-paragraph warning** that
   the EMC baseline constant is quoted against CT18ANLO and that any future
   EMC hook must state its own proton denominator, so the two do not silently
   drift. That is cheaper now than a fix round later.

---

## 7. Adjacent defects found while mapping (not part of D2's brief)

Three, all verified against the build. They are reported here rather than
fixed, because D2 is read-only and none of them is the injection map.

**(i) `PipelineConfig::kernel` is silently ignored on tagged channels, and the
`meta` then lies about it.** `Pipeline` reads `cfg_.kernel` only in the
non-tagged branch (`src/core/pipeline.cpp:951` (as of adec442)); the tagged branch
(`:853-854`) goes straight to `make_struck_cluster_source`. `validate()` refuses
`kernel` only together with a non-Miller `b1_model` (`:555-560`), so a tagged
config with a kernel attached passes validation. Measured:

```
tagged-6Li-alpha run with cfg.kernel = InclusiveKernel(d, f2_source=CT18NLO):
  validate() raised nothing
  the sampler actually used  F2A(0.3, 10) = 0.36911493451624272   <- ToyF2
  CT18NLO NuclearF2(d) would be              0.46397034417540034
  meta['b1_model'] = 'caller-supplied kernel'
  meta['b1_unpol'] = 'caller-supplied kernel'
```

So the npz records a provenance for a kernel that never ran, and the numbers are
the toy ones. This is the exact failure mode the b₁ provenance rules exist to
prevent, in a corner they do not cover. **Fix: `validate()` should refuse
`kernel` on any non-inclusive channel** (it reaches neither the tagged nor the
coherent rate), or `Pipeline` should honour it on the tagged branch. Refusing is
simpler and matches the existing rule. Note that closing §5.8 clause 3 does
*not* close this one — this is `kernel` alone, with no new flag set.

**(ii) An `--isotope d --channel inclusive` run gets the ⁶Li rank-2 transfer.**
`default_inclusive_kernel`'s tensor branch keys on `ion.spin == 1`
(`src/core/pipeline.cpp:140`), and the deuteron is spin 1, so it receives
`Li6B1(MillerB1)` — i.e. `LI6_B1_RANK2_TRANSFER × LI6_B1_PER_NUCLEON` =
0.921947 × 1/3 = **0.307316** applied to a deuteron's own b₁. Measured:
`default_inclusive_kernel(deuteron()).tables(x, q2).b1` is *identical* to the
⁶Li kernel's at every point, and is 0.307316 × `toy_b1(x, q2, f1)`. The isotope
guard exists in `validate()` (`:571-589`) but only inside the **non-Miller**
branch, and Miller is the default. `docs/PHYSICS_CHANNELS.md:128` already
describes this behaviour ("for ANY spin-1 ion"), so it is documented — but a
deuteron run is not a ⁶Li run and the 2/6 embedded-deuteron counting factor has
no meaning for A = 2. Worth a decision in phase F.

**(iii) `PythiaBridgeOptions::f2_source` is unbound.** §5.6. Harmless today,
a live inconsistency the moment `--unpol-sf` lands.

---

## 8. Effort, honestly

**What was done:** the tree at `903fcc9` was read end to end for every
construction and defaulting site of `InclusiveKernel`, `UnpolSF`, `PolSF`,
`TensorSF` and `RFunc` across `src/core`, `src/pythia`, `src/lhapdf`,
`include/`, `examples/`, `python/bindings.cpp`, `python/lipolgen/` and
`validation/`; **every** `file:line` in this document was re-read after the
text was written — 95 citations, 77 distinct, all of which resolve to an
existing in-range line, and each of which was printed and checked to land on
the construct the surrounding sentence claims. The
Phase-A `B1UnpolSource` thread was read in full across
`pipeline.hpp` / `pipeline.cpp` / `bindings.cpp` / `__init__.py` / `cli.py` /
`test_b1_model.py`. The installed LHAPDF store was listed and the three sets'
`.info` metadata read.

**What was measured, on this build, today:** §3a (four sampler builds, two
targets), §3b (two kernel builds × 7 points × 2 observables, plus the g₁ⁿ
sign table), §3c (two kernel builds × 4 points), §3d (3051-cell R statistics),
§3e (grid limits and below-grid F₂ probes), §7(i) (one tagged pipeline run with
a kernel attached, plus its npz meta), §7(ii) (two `default_inclusive_kernel`
calls). Every number in this document that is presented as a measurement was
produced by one of those runs; every number presented as a quotation from the
tree is a quotation.

**What was NOT done, and is therefore not claimed:**

* **No implementation.** No enum, no field, no binding, no flag, no test was
  written. §5 is a design, not a diff, and its line counts and refusal messages
  are proposals.
* **The `--pol-sf` measurements are `LhapdfG1` construction only.** They were
  taken by hand-building `InclusiveKernel::Options::g1_model`, which is the
  object the design would install, but **no pipeline was run end to end on a
  polarised backend** and no asymmetry was extracted. The A_∥ ratios in §3b are
  kernel-level, at one y, not run-level.
* **The tagged rate ratio in §3a is the isoscalar deuteron only.** The triton
  (⁷Li α-tag) and `NEUTRON_TARGET` (d+p tag) `dis_target`s were not measured,
  and they are exactly the two where the toy n/p ratio and the toy g₁ⁿ sign
  (§3b) do the most damage. Somebody should measure them before the flag is
  quoted on those channels.
* **MSTW's grid boundary was not characterised** the way CT18NLO's was; the
  42.33 % below-grid figure is CT18NLO's alone. MSTW's below-grid F₂ᵖ was
  probed at two points and looked mild, which is an observation and not a
  measurement of the fraction.
* **Nothing was measured about EPPS21 as a *kernel* backend**, because no
  kernel can take it: `Epps21Ratio` is not an `UnpolSF` and there is no
  `emc_ratio` wiring to test.
* **The `--unpol-sf`/`--pol-sf` orthogonality claim of §3c was checked at four
  points at one Q² each.** It is a demonstration that A₁ moves, not a
  characterisation of by how much across the window.
