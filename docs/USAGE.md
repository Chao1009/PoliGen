# LiPolGen — the C++ API in two pages

Everything below is `namespace lipolgen`, C++17, header `<lipolgen/pipeline.hpp>`,
link `LiPolGenCore` (+ `LiPolGenHepMC` if you write HepMC3).

```bash
source env.sh
cmake -S . -B build && cmake --build build -j8 && ./build/lipolgen_tests
```

### Installing the Python package with pip

The in-tree build above is what `env.sh` and every example on this page
assume. If you just want `import lipolgen` in some other virtualenv, without
sourcing `env.sh` or touching `build/`, use the `pyproject.toml`
(scikit-build-core) instead:

```bash
LIPOLGEN_DEPS_PREFIX=/path/to/deps/install pip install -e /path/to/LiPolGen
```

`LIPOLGEN_DEPS_PREFIX` (HepMC3/LHAPDF/PYTHIA8, same meaning as in `env.sh`)
can also be passed as `--config-settings=cmake.define.LIPOLGEN_DEPS_PREFIX=...`;
the environment variable is a convenience CMakeLists.txt reads itself when
the CMake cache variable is not already set. Two of the three dependencies
still need their own data at run time even once the module is installed —
export these the same way `env.sh` does, since `pip install` does not:

```bash
export PYTHIA8DATA=$LIPOLGEN_DEPS_PREFIX/share/Pythia8/xmldoc
export LHAPDF_DATA_PATH=$LIPOLGEN_DEPS_PREFIX/share/LHAPDF
```

(LiPolGen's own `data/vmc` tables are small enough that the wheel vendors
them directly, so no third variable is needed for those — see
`data_dir()`/`$LIPOLGEN_DATA_DIR` in `include/lipolgen/cluster.hpp`.)

**Portability caveat.** The extension's RPATH is baked in at build time as
`$ORIGIN` (LiPolGen's own libraries, shipped alongside it in the wheel) plus
the literal `LIPOLGEN_DEPS_PREFIX` path (HepMC3/PYTHIA8/LHAPDF, NOT shipped
in the wheel) — so a wheel built this way only runs on the machine it was
built on, at that same path.
[`auditwheel repair`](https://github.com/pypa/auditwheel) fixes the *library*
half of that, and since 2026-09-05 that is measured rather than recommended:
run on this tree it vendors `libHepMC3.so.4`, `libLHAPDF.so` and
`libpythia8.so` into a `lipolgen.libs/` directory, rewrites the RPATHs to
`$ORIGIN`-relative, and takes the wheel from **1 779 767 to 7 232 971 bytes**
under the tag **`manylinux_2_35_x86_64`**. The repaired wheel was verified to
import, generate events and run the T2 tier with the dependency prefix
*removed from the filesystem altogether*.

It does **not** fix the *data* half. PYTHIA's `xmldoc` and LHAPDF's PDF sets
are not in the wheel and their compiled-in defaults are the build machine's
absolute paths, so on any other machine `PYTHIA8DATA`, `LHAPDF_DATA_PATH`
(and `LIPOLGEN_PYTHIA8_PDFDATA`, for `MstwSF`) must be exported exactly as
above — with them unset the failures are
`Couldn't find required lhapdf.conf system config file` and
`PYTHIA Error in Settings::init: settings file <build machine path>/Index.xml
not found`. Every command, every number and every one of those messages is in
**`docs/PACKAGING.md`**.

Repairing the wheel also turns it into a genuine combined/linked work with all
three (GPL-2-or-later, GPL-3.0, GPL-3.0), so whoever redistributes an
`auditwheel`-repaired wheel is redistributing under GPL-3.0-or-later terms for
the combination, regardless of LiPolGen's own license header (see
`docs/OPEN_ITEMS_SOLUTIONS.md` §12–13 and `docs/open_items/engineering.md` §C).

## 1. The shape of a run

A run is **one `PipelineConfig` + one `RunPlan`**.

* `PipelineConfig` says *what physics and how much of it*: isotope, beam
  configuration index (0/1/2 = low/mid/top), channel, acceptance window
  (`Scenario`), far-forward optics, seed, and either `lumi_pb` or `n_events`.
* `RunPlan` (`bookkeeping.hpp`) says *how the luminosity divides between spin
  fills*: `tensor_thirds_plan`, `helicity_flip_plan`, `transverse_tensor_plan`,
  `tensor_flip_plan`, or your own `std::vector<SpinCategory>`.

**Fixed-count mode or luminosity mode, and how the CLI reaches each.**
`--events N` fixes the total; `--lumi X` [pb⁻¹] derives the counts from the
luminosity and the cross sections. They are exclusive, and **`--lumi X` alone
is enough** — since 2026-09-05 typing it with no `--events` anywhere sets
`events = 0`, which is what selects luminosity mode. (Before that the option
defaults carried `events = 100000`, so `--lumi X` alone was refused with
"`--events and --lumi are exclusive`" and the mode was reachable only as
`--events 0 --lumi X`, which no help text said. An `--events` typed on the
command line or given in a `--config-file` still collides with `--lumi` and is
still refused.)

`Pipeline` resolves both in its constructor — grids, per-spin-state amplitude
tables, spectator amplitude tables, per-category counts — and then maps
**event index → `Event`** as a pure function. That is the whole design:

```cpp
Event Pipeline::event(std::uint64_t index) const;   // pure, thread-safe
```

Category `k` and local counter `j` are recovered from the index and the stream
is `Rng(seed, run, k, j)`. Nothing is carried between events, so
`for_each(sink, nthreads)` is bit-identical to `for_each(sink)` for any thread
count, and `event(i)` can be called from anywhere, in any order.

Three ways to consume a run:

```cpp
p.for_each([](const Event& ev) { ... });          // streaming, nothing stored
p.for_each(sink, /*nthreads=*/8);                 // same events, same order
p.for_each_range(sink, 1000, 2000);               // streaming over a range
Event ev;  while (p.next(ev)) { ... }             // pull API (cursor; rewind())
std::vector<Event> batch;                          // explicit index range
p.generate_range(1000, 2000, batch);              // resized, reconstructed in place
```

`event(index, ev)`, and therefore all of the above, reconstructs into `ev`
**in place** (`Event::reset()` keeps the record's heap capacity), so a loop
over a reused `Event` allocates nothing per event.

Bookkeeping, all share-invariant in pb:

```cpp
p.sigma_per_category_pb();   // [pb], no lumi_fraction, no Scenario::run_share
p.sigma_pb();                // sum_k lumi_fraction_k * sigma_k
p.lumi_per_category_pb();    // [pb^-1] (luminosity mode), x Optics::lumi_fraction
p.counts();  p.size();       // events per category / in total
p.optics_lumi_factor();      // the Optics::lumi_fraction that reached the counts
```

**The optics luminosity fraction reaches the COUNTS, never the cross
sections.** A far-forward working point that buys acceptance by de-squeezing
β*_x pays for it in luminosity: `Optics::lumi_fraction` is 1 at the Yellow
Report envelopes and 0.1467 at the ⁶Li 5×41 tagging point. It multiplies the
per-category luminosity (`PipelineConfig::apply_optics_lumi_fraction`, default
`true`), so at a fixed `lumi_pb` the tagging optics deliver 0.1467× the events
of the YR optics — 1 116 408 against 7 612 479 in the gate — while
`sigma_per_category_pb()` is bit-identical between the two. Set the flag to
`false` to quote a yield at the machine luminosity regardless of the optics.

## 2. Inclusive

```cpp
#include "lipolgen/pipeline.hpp"
using namespace lipolgen;

PipelineConfig cfg;                       // defaults: 6Li, config 1, generator window
cfg.channel  = PipelineChannel::Inclusive;
cfg.isotope  = "6Li";
cfg.lumi_pb  = 2.0;                       // or cfg.n_events = 1'000'000;
cfg.seed     = 20260713;
// cfg.kernel = my_kernel;                // default: Miller b1 through the 6Li
                                          // rank-2 transfer + 1e-2 toy Delta,
                                          // and target_mass ON (xsec.py's own
                                          // default since 2026-08-29)
Pipeline p(cfg, tensor_thirds_plan(0.7, 0.6));

p.for_each([](const Event& ev) {
  const Particle* e1 = ev.find(Role::ScatteredElectron);
  printf("x=%g Q2=%g phi=%g  m=%g  E'=%g\n",
         ev.kin.x, ev.kin.q2, ev.kin.phi, ev.spin.m_ion, e1->p.e);
});
```

Record: beam e⁻, beam ion, e′, struck nucleon (status 3), γ\* (status 3),
X (pdg 92). The balance is **per-nucleon**, `k + P_N = k' + X`, exactly as
`generator.hpp` defines it — the (A−1) remnant is not written, so the
whole-nucleus balance and the total charge are deliberately open here. This path
is bit-identical to `InclusiveGenerator::run_n` for the same seed.

The struck nucleon is on shell at the **free** nucleon mass `M_NUCLEON` and its
species is drawn `Z F2p(x, Q²) : N F2n(x, Q²)` on the kernel's own unpolarized
backend — not flat `Z : N`, which is the x-independent limit of it (0.5 against
a true 0.6154 for ⁶Li at x = 0.5). `GeneratorConfig::struck_nucleon_pdg` pins
the species; `InclusiveGenerator::proton_fraction(x, q2)` is the rule itself.

## 2a. b₁ backends for the inclusive channel — `--b1-model`

The inclusive kernel's rank-2 slot is filled by `default_inclusive_kernel`
whenever `PipelineConfig::kernel` is null, and **which** b₁ goes in there is a
choice, not a fact:

| `--b1-model` | `B1Model` | what it is |
|---|---|---|
| `miller` (default) | `B1Model::Miller` | `Li6B1(MillerB1)` through `LI6_B1_RANK2_TRANSFER` and `LI6_B1_PER_NUCLEON` (2/6). **Bit-for-bit what every published inclusive tensor number was made with**, and pinned at rtol 1e-12 by `validation/reference/b1_default_li6.json` (`tests/test_b1_nuclear.cpp` T9). This IS the run PLAN's "toy": today's default `b1_func` is `Li6B1(MillerB1)`, which is what `toy_b1` reaches. |
| `cdks` | `B1Model::Cdks` | the same ⁶Li rank-2 transfer on the **other camp** for b₁ᵈ — the digitized Cosyn–Dong–Kumano–Sargsian PRD **95** (2017) 074036 Fig. 4 column. \|b₁\| is **two orders of magnitude smaller below x ≈ 0.1 and comparable to or larger than Miller's above it**: peak \|x·b₁\| **3.34e−4** against Miller's 4.27e−4 over x ∈ [0.02, 0.95] at Q² = 2.5, and at x = 0.3 it is **8× larger** with the opposite sign (∫b₁ dx **1.39e−4** against 9.1e−4). Different sign structure throughout. **These numbers DOUBLED on 2026-09-03**: `b1_convolution()` stopped applying `B1_PER_DEUTERON_TO_PER_NUCLEON` to a column that CDKS Eq. (10) and the text under their Eq. (16) say is already per nucleon (§2a "The two b₁ camps are not per the same thing"). |
| `li6-convolution` | `B1Model::Li6Convolution` | `Li6ConvolutionB1` (`b1_nuclear.hpp`): the **four-term α–d convolution** of `docs/open_items/run_2026-09-02/design_D_b1_li6.md`. |

Miller (HERMES-like) and CDKS (convolution) are *different camps* for the
deuteron's own b₁ and the library does not adjudicate between them. **Say which
one a plot used.**

### The two b₁ camps are not per the same thing

Both digitized curves are consumed as **per-nucleon** b₁, because every F₁ in
this program is per nucleon — but they do not arrive that way, and since
2026-09-03 they no longer share one conversion constant
(`include/lipolgen/constants.hpp`).

| camp | table | applied factor | why |
|---|---|---|---|
| Miller | `kB1Miller` | `B1_MILLER_TABLE_TO_PER_NUCLEON` = **0.5** | Miller PRC **89** (2014) 045203 is **per deuteron**: his Eq. (1) densities are "in a target hadron", Eq. (5) is a light-cone correlator in the normalised deuteron state, and Eq. (6)'s ½ is the quark-spin average of a *spinless* pion, so no 1/A is left anywhere in Eqs. (1)/(5)/(6)/(20). **Likely, not certain** — see below |
| CDKS | `kB1CdksQ2p5` | `B1_CDKS_TABLE_TO_PER_NUCLEON` = **1** | CDKS PRD **95** (2017) 074036 is **already per nucleon**: Eq. (10)'s spectral function carries an explicit 1/A, the text under Eq. (16) says so in words, f(y) is normalised to one nucleon and F₁ᴺ = (F₁ᵖ + F₁ⁿ)/2. **Certain** |

The arbiter is neither paper but **HERMES**, the data both camps plot against:
their published b₁ᵈ is **per nucleon**, because their Eq. (5) divides by an F₁ᵈ
built from F₂ᵈ = (F₂ᵖ + F₂ⁿ)/2 — and inverting their own Table II reproduces
that F₁ᵈ in all six bins (mean ratio 0.95, against 0.47 for the per-deuteron
one). Full argument:
`docs/open_items/run_2026-09-03/phase_A_miller_normalisation.md`.

**What changed at the run surface:** `--b1-model cdks` (and anything reading
`b1_convolution` / `CdksB1`) is a **factor 2 larger** than before 2026-09-03.
`--b1-model miller`, the default, and `--b1-model li6-convolution`, which
already reached the raw column through `cdks_b1_raw_per_nucleon()`, are both
**bit-for-bit unchanged**; the rtol-1e−12 reference gate does not move.

**What is still open on Miller.** His paper is self-inconsistent by exactly
this factor: the derivation is per deuteron, but his Table I transcribes
HERMES's per-nucleon numbers unrescaled, his Fig. 5 overlays them on the curve,
and he tunes P₆q to one of them. Keeping the 0.5 is the status quo and is
recorded as an author decision (`OPEN_ITEMS_SOLUTIONS.md` §10, condition 6).
The sharp form of the question is *"does Miller's Eq. (20) evaluate to
10.5 × 10⁻² or 5.25 × 10⁻² at x = 0.012?"*, and answering it needs the pion PDF
set (his ref [29], model 1) that this tree does not have.

```cpp
PipelineConfig cfg;
cfg.channel  = PipelineChannel::Inclusive;
cfg.isotope  = "6Li";                      // 6Li ONLY for li6-convolution
cfg.n_events = 1000000;
cfg.b1_model = B1Model::Li6Convolution;
cfg.b1_band_scale = 1.0;                   // run 0, 1 and 2 -- see below
cfg.b1_alpha_d_dwave_weight = 1.0;         // terms (2d)+(2a) knob
cfg.scenario.x_max = 0.95;                 // needed; see "the top x cell"
```

```
lipolgen-run --channel inclusive --isotope 6Li --events 1000000 \
             --b1-model li6-convolution --b1-band-scale 1 --x-max 0.95
```

Every generated file records what made it: `meta["b1_model"]`,
`meta["b1_band_scale"]` and `meta["b1_alpha_d_dwave_weight"]` are in the npz /
HFS metadata unconditionally, so three otherwise identical runs are not
indistinguishable. A caller-supplied `cfg.kernel` is labelled
`"caller-supplied kernel"` rather than with a flag that did not run.

### Where the model is legal, and why

`validate()` refuses **either opt-in model** — `cdks` as well as
`li6-convolution`, because `cdks` is `Li6B1`'s ⁶Li rank-2 transfer too —
outside **the inclusive channel** and outside **isotope `6Li`**, each with its
own message. Elsewhere `default_inclusive_kernel` would never read the flag,
while the npz/HFS `meta` would still record it.

* **Inclusive only.** `Li6ConvolutionB1` models the b₁ of the whole ⁶Li *as
  seen inclusively*: the z-smearing of the embedded-deuteron term and the
  orbital alignment of the α–d D-wave terms are both properties of the α–d
  relative motion. On a **tagged** channel `TaggedSampler` already draws the
  α–d momentum *and* its m-dependent angular correlation from n_M(k, k̂), so
  that same wave function is in the event weight already — this is the double
  counting `StruckClusterOptions::inclusive_b1` warns about, and folding the
  convolution in on top would count it **three** times. On the **coherent**
  channel the tensor signal lives in the recoil azimuth and the inclusive b₁ is
  not folded in at all.
* **⁶Li only, keyed on the ISOTOPE and not on the spin.** The inclusive channel
  accepts `isotope = "d"`, which is spin 1, so a spin test would let a deuteron
  beam silently run the ⁶Li α–d convolution — N_αd, the α–d densities, the 2/6
  and 4/6 counting factors — on a deuteron. For the A = 2 kernel use
  `DeuteronConvolutionB1` directly; that is the validation gate, not a beam
  species. `7Li` is spin 3/2 and has **no rank-2 input here — see §2c, which
  is about running ⁷Li and not about this flag**.
* A **caller-supplied `cfg.kernel`** together with a non-`miller` `b1_model` is
  also refused: the kernel wins and the flag would be silently ignored, so the
  two are saying contradictory things.
* **A knob that the chosen backend does not read is refused too.**
  `--b1-band-scale` is deliberately not applied to `miller` (those numbers are
  the published ones) and only `li6-convolution` reads
  `--b1-alpha-d-dwave-weight`, but all three keys go into the file's metadata
  unconditionally — so `--b1-model miller --b1-band-scale 0` or
  `--b1-model cdks --b1-alpha-d-dwave-weight 2` would record a variation that
  never ran, and both now throw. Both knobs must also be ≥ 0 on every model.

### Which R each object uses — they are NOT the same, and that is decided

`Li6ConvolutionOptions::r_func` null means **`r_sigma_lt`**, the kernel's own
toy R, so that the ⁶Li backend's F₁ is consistent with `InclusiveKernel`'s.
`DeuteronConvolutionB1::Options::r_func` null means **`r1998`**, the SLAC world
fit, because that is what CDKS used and the A = 2 gate exists to reproduce
*their* figure. So **the R the gate was validated with is not the R the shipped
⁶Li backend runs with.**

**Both defaults stay, and it is an author decision rather than an oversight**
(`OPEN_ITEMS_SOLUTIONS.md` §10, condition 4). What decides it is that the
observable is a *ratio*: the tensor weight is K/D_φ with K ∋ b₁ and D_φ ∋ F₁,
and the F₁ in the denominator is `InclusiveKernel`'s, built from the shared
`ToyF2` with **its** R — which is `r_sigma_lt`, because nothing sets it. Give
`Li6ConvolutionOptions` a different R and the (1 + R) in the numerator no
longer cancels the one in the denominator: at Q² = 2.5 the mismatch factor
(1 + r1998)/(1 + r_sigma_lt) is 1.088 at x = 0.1 and 1.039 at x = 0.3 — the
same size as the whole effect being chased. Fidelity to CDKS is worth having
on the *gate*, which has no denominator; inside the pipeline consistency wins.

**The measured cost of the choice** (Q² = 2.5, default options, 2026-09-03):

| x | x·b₁ with `r_sigma_lt` (the default) | with `r1998` | change |
|---|---|---|---|
| 0.05 | +3.650716e−6 | +3.449387e−6 | −5.5 % |
| 0.10 | −2.796612e−6 | −3.484857e−6 | +24.6 % |
| 0.20 | −2.120522e−5 | −2.235704e−5 | +5.4 % |
| 0.30 | −4.173248e−5 | −4.284379e−5 | +2.7 % |
| 0.50 | +5.173934e−5 | +5.075777e−5 | −1.9 % |

**Where that comes from, and it is not the (1 + R) prefactor.** Term (1), the
embedded deuteron, is **bit-identical** under the swap — it carries b₁ᵈ from
the injected `TensorSF` (the raw digitized column), which has no R at all. The
whole effect is in the two orbital terms, whose F₁ᵈ slot *is* `f1_cdks`, and
there it is **×0.54 to ×0.94**, far more than the ~8 % the prefactor allows:
those terms convolve F₁ᵈ(x/z) against a density that **integrates to zero**, so
they respond to the *slope* of R, and `r_sigma_lt` is x-independent by
construction while `r1998` runs from 0.30 at x = 0.05 to 0.20 at x = 0.5. The
+24.6 % at x = 0.10 is that ×0.64 seen through the near-cancellation between
term (1) and terms (2d)+(2α). On the **gate**, where there is no cancellation,
the same swap is worth only +2.7 % on G3b (0.8432 → 0.8662 on a 0.01 grid over
[0.10, 0.80], MSTW at Eq. (21)) — so the gate's verdict does not rest on it.

**The third option now EXISTS, opt-in — `--r-source` (2026-09-06).**
`--r-source {unset, sigma-lt, r1998}` (`PipelineConfig::r_source`, `RSource`)
threads **one** R object through `default_inclusive_kernel` into both
`Li6ConvolutionOptions::r_func` and `InclusiveKernel::Options::r_func`, so the
ratio's two halves cannot disagree by accident. `unset` is the default and
**does not enter that branch at all** — today bit for bit, by construction —
and the flag is refused off `--b1-model li6-convolution`, where the numerator
has no F₁ of its own to give an R to. `sigma-lt` installs the shared object
and is measured bit-identical to `unset`, which is the wiring's own test and
is why its provenance row reads `not-read`.

**It is a wiring change, but its cost is a physics number, and it does not
simply undo the mismatch** (`docs/open_items/run_2026-09-06/phase_A_numbers.md`
§A1, which carries the full six-point × three-y table):

| at y = 0.5, Q² = 2.5 | x = 0.05 | x = 0.10 | x = 0.30 |
|---|---|---|---|
| K/D_φ (= −3/2 · A_zz), `r1998` in the **numerator only** | −5.5148 % | +24.6099 % | +2.6629 % |
| K/D_φ, `--r-source r1998` (**both halves**) | −3.8357 % | +26.3988 % | +3.3518 % |
| cos 2φ amplitude, numerator only | 0 (bit-identical) | 0 | 0 |
| cos 2φ amplitude, `--r-source r1998` | −8.3268 % | −6.7268 % | −3.1441 % |

Three things to read off it. (a) Sharing the hook makes the shift **larger**
at x = 0.10, so it restores part of the cancellation and not all of it — the
orbital terms are not proportional to F₁. (b) The cos 2φ amplitude separates
the two options **exactly**: Δ carries no R, so a numerator-only swap cannot
move it and the shared hook moves it purely through D_φ. (c) The shared
shift is **y-dependent** where the numerator-only one is not (it agrees to
2 ulp across y): at x = 0.05, Q² = 2.5 it runs −5.4705 % at y = 0.1 to
+2.3678 % at y = 0.9, **straddling zero** — so a `--r-source r1998` number
must be quoted with its y. And unlike `--b1-unpol`, this flag moves the
**unpolarised rate** too (σ −0.684681 % on a 2000-event ⁶Li run at
`--x-max 0.95`), because the kernel's `r_func` reaches F₁, F_L, D(y) and
`ToyG1`.

**The default is unchanged and nothing here is decided**: registry row 3 is
still open, and the option it called "the better third one" is now priced
rather than taken.

### The mandatory 100 % band — never quote a single row

`--b1-band-scale {0,1,2}` multiplies the **whole** b₁ (all four terms of
`li6-convolution`, and the whole `cdks` b₁; `miller` is deliberately untouched,
because those numbers are the published ones). The rationale is a measurement:

> **Q(⁶Li) = −0.0818(17) fm²** against **Q_d = +0.2859(3) fm²** — the α–d
> relative D wave enters the closest measured observable with the **opposite
> sign** to the deuteron's own D state and very nearly cancels it.
> *(Q(⁶Li) is the repository's single copy, `LI6_QUADRUPOLE_FM2` in
> `include/lipolgen/rc.hpp`: TUNL's A = 6 evaluation, 1998CE04. Pyykkö's
> compilation, Mol. Phys. **106** (2008) 1965, gives −0.0806(6) fm² from the
> molecular-beam measurement of Cederberg et al., Phys. Rev. A **57** (1998)
> 2539 — 1.5 % away, and nothing here depends on the choice. Q_d: Bishop &
> Cheung, Phys. Rev. A **20** (1979) 381.)*

The two are different operators (charge quadrupole vs light-cone momentum
alignment) so the cancellation need not carry over, and the sign flip **is real
in the VMC overlaps** — but it is k-dependent: φ₀φ₂ is opposite to the
deuteron's below the S node at 0.678 fm⁻¹ (66 % of the density) and
deuteron-like between the nodes (33 %), so the net sign of the orbital terms is
a computed output, not an assumption. Therefore: **every published number from
these backends is {0, 1, 2} × b₁.** `--b1-alpha-d-dwave-weight 0/1/2` is the
*shape* variant of the same worry and is reported separately.

### The A = 2 validation gate — PASSED since 2026-09-03, and what that does not say

Design §5's gate is **blocking**: fed the AV18 deuteron u(k), w(k) instead of
the α–d waves, the same kernel must reproduce the digitized CDKS Fig. 4 before
any ⁶Li number is quoted. **It now does.** The verdict is read where design
§5.4's checklist ends — with **MSTW2008 LO, CDKS's own nucleon PDF**
(checklist item 4; `MstwSF`, below), at **CDKS Eq. (21)'s δ-function**, which
is `DeuteronConvolutionB1`'s default:

| clause | result |
|---|---|
| G0a–G0g (analytic identities), G1a–G1f (unpolarized convolution), layer 2 (α–d sign gate) | **PASS** |
| **G3a** — exactly two sign changes in (0, 1.0], the first falling within ±0.08 of 0.0656, the second rising within ±0.10 of 0.4572, the peak within ±0.10 of 0.766 | **PASS, QUALIFIED**: zeros **0.0279** and **0.4952**, peak at **0.7716** — misses of 0.038 / 0.038 / 0.0059 against tolerances 0.08 / 0.10 / 0.10. The qualification is the window's **upper** edge, which no tolerance derives: over the reference's own domain [0.010, 1.590] the computed curve has a **third** sign change (1.2204 toy / 1.2177 MSTW / 1.1977 CT18NLO) and the reference has none there. Never quote "exactly two" without the window; the numbers and the argument are in `docs/OPEN_ITEMS_SOLUTIONS.md` §10, "G3a's stated limitation" |
| **G3b** — max \|x·b₁\| over [0.10, 0.80] within a factor 2 of the digitized 1.08521e−3 | **PASS**: **9.15096e−4**, ratio **0.843243**, a factor **1.19 low** — inside [0.5, 2] with the lower edge cleared by a factor 1.69 |
| G3c — Close–Kumano ∫b₁ dx, **reported, not enforced** | +2.24896e−4 against the digitized +4.59200e−4, i.e. 0.49 of it (with CD-Bonn, +4.48580e−4 = 0.98) |

**Four things that belong in the same breath as the pass.**

1. **It is conditional on the CONFIGURATION, not on the build.** The gate
   passes for the **MSTW2008 LO** nucleon input at CDKS Eq. (21)'s δ-function.
   The **shipped default** unpolarised backend — the library `ToyF2` — gives
   ratio **0.440** on the same clause, *outside* the [0.5, 2] window. MSTW is
   checklist item 4 and the dominant term of the design's own residual budget,
   so reading G3b off `ToyF2` reads it *before* the checklist. Since
   2026-09-04 the passing configuration is emittable from the run surface:
   `--b1-unpol mstw` (below) puts `MstwSF` into
   `Li6ConvolutionOptions::unpol` through the pipeline's own
   `default_inclusive_kernel`. So the rule is **quote numbers made with
   `--b1-unpol mstw`; the default `toy` backend sits outside the gate's
   acceptance window and its numbers are not covered by the lift.** They
   cannot be rescaled into covered ones either — the gap is a shape change,
   mstw/toy = **1.848 / 1.276 / 0.817** at x = 0.10 / 0.30 / 0.50, up to a
   factor 1.85 and not monotone. Selecting `mstw` needs the optional PYTHIA
   tier and is **refused at configuration time, never downgraded**, when the
   tier or the `mstw2008lo.00.dat` grid is missing; in such a build the
   doctest's MSTW rows are skipped **loudly**, and since 2026-09-05 that is
   true of every row the lift rests on: the verdict row `T1v`, the
   gate-condition-3 row `item 5v` (CD-Bonn + MSTW, 1.000338) and the four
   `MstwSF` backend cases carry `doctest::skip(!mstw_grid_present())`, so
   **doctest reports them as skipped, not as passed**. Two rows that are not
   themselves gate conditions — T16's `--b1-unpol mstw` reach row and T17's
   third-crossing row, each inside a case whose shipped-default rows must
   still run — keep an in-case skip: they print
   `SKIPPED (MSTW2008 LO unavailable)` with what was not measured and are
   tallied as passed (the choice between splitting them out and leaving them
   is `AUTHOR_DECISIONS.md` B26). **They are no longer the only two rows in
   the tree with an in-case skip**, but they are the only ones left that are
   *tallied as passed*: the twelve `validation/reference` port gates in
   `test_reference.cpp`, `test_bookkeeping.cpp`, `test_spectator.cpp` and
   `test_tagged.cpp` kept the same `MESSAGE`+`return` shape — and so passed
   with **zero assertions** on a checkout without the blobs — until
   2026-09-16, when they were given the same registration-time
   `doctest::skip` and a `REQUIRE` on the loader (a blob that IS there and
   does not parse is now a failure, not a skip). Either way a build that
   cannot reproduce the verdict row also cannot emit a number that claims it.
2. **It passes comfortably at Eq. (21) and marginally at Eq. (17).** At κ = 1
   MSTW gives **0.520** — inside [0.5, 2] by 4 % of its own value. Any
   statement of the form "the gate passes" that does not also say "at CDKS
   Eq. (21)'s δ-function" is overselling it.
3. **With CD-Bonn as well the residual closes.** `Options::wave = kCdBonn`
   (below) — CDKS's own wave function — gives peak ratio **1.000338**, zeros
   0.0641 / 0.4570 against the digitized 0.0656 / 0.4572, and a dip of
   −1.769065e−4 against −1.768140e−4. Read that as *"the residual is now below
   the error of digitizing a published figure"*, never as three-digit agreement
   with CDKS; and it is specific to CD-Bonn **and** MSTW together — on CT18NLO
   the same swap improves two landmarks and degrades two others.
4. **The gate is A = 2.** It validates the kernel on the **deuteron**. There is
   no measured b₁ for any A > 2, so **the mandatory ±100 % band on every ⁶Li
   number stays** — that rule comes from Q(⁶Li) = −0.0818(17) fm² against
   Q_d = +0.2859(3) fm², not from the gate.

**Since 2026-09-03 the gate is quoted at CDKS Eq. (21)'s δ-function**
(`DeuteronConvolutionB1::Options::finite_q_delta`, `true`): Eq. (21) is their
exact definition, y = (E − p_z κ)/M_N with κ = \|q⃗\|/ν, and Eq. (17)'s
(E − p_z)/M_N is the "≃" of their Eq. (18). It is worth ×1.62 at the peak.
`Li6ConvolutionOptions::finite_q_delta` keeps its `false` default, where the
same switch is worth 1 %.

**The residual budget, all measured** (doctest grid, `linspace(0.001, 1.59,
300)`), starting from the κ = 1 `ToyF2` column:

| item | effect on the peak ratio | in the default? |
|---|---|---|
| starting point: `ToyF2`, κ = 1 | 0.271689 | — |
| 0. finite-\|q⃗\| δ-function, κ = 1 → √(1+γ²) | ×1.6195 | **yes** |
| 4. nucleon PDF, `ToyF2` → **MSTW2008 LO** | ×1.9152 | opt-in (`MstwSF`) |
| 5. wave function, AV18 → **CD-Bonn** | ×1.1863 | opt-in (`wave = kCdBonn`) |
| 1. target mass (CDKS Eq. 22) / 2. R (`r1998`) | ×1.49 at x = 0.8 / ×0.98 | yes / yes |
| **AV18 + MSTW + Eq. (21)** | **0.843243** | the verdict row |
| **CD-Bonn + MSTW + Eq. (21)** | **1.000338** | — |

Items 0 and 4 are independent to 0.07 % (0.271689 × 1.91519 × 1.61951 =
0.842682 predicted against 0.843243 measured), so there is no third
unexplained factor at the peak, and item 5 supplies exactly the 1.186 that
was left over. **The ⁶Li publication ban that the previous failure imposed
under design §5.4's *Escalation* clause is therefore lifted** — that clause
reads "if **after the checklist** the peak ratio is still outside a factor of
2", and after the checklist it is not. Quote ⁶Li numbers **as a band**, and say
which configuration made them.

Full measurements: `docs/open_items/run_2026-09-03/phase_A_numbers.md` (the
rerun and the CD-Bonn section) and `phase_A_cdbonn.md` (the coefficients). The
previous, failing verdict is `docs/open_items/run_2026-09-02/phase_D_gate.md`,
**superseded on G3b** and kept for the argument, not for the number.

#### `MstwSF` — MSTW2008 LO from the grid PYTHIA already ships

`include/lipolgen/mstw_sf.hpp`, implemented in `src/pythia/mstw_sf.cpp`; the
core library still links neither LHAPDF nor PYTHIA. CDKS computed their b₁ with
MSTW2008 LO, and the CENTRAL member is on disk as
`<pythia8 datadir>/pdfdata/mstw2008lo.00.dat` — `Pythia8::MSTWpdf` with
`i_fit = 3` reads exactly it, so the long-standing "MSTW2008 LO is not
installed" in the older documents was true only of LHAPDF's set store.

```cpp
#include "lipolgen/mstw_sf.hpp"
DeuteronConvolutionB1::Options o;
o.unpol = std::make_shared<const MstwSF>();   // i_fit = 3 = MSTW2008 LO
DeuteronConvolutionB1 gate(o);
```

```python
o = lipolgen._lipolgen.DeuteronConvolutionB1.Options()
o.unpol = lipolgen._lipolgen.MstwSF()
```

The charge weights and the F₂ construction are copied verbatim from
`LhapdfSF`, so an MSTW/CT18NLO ratio taken from these two classes is a **PDF**
comparison and not a convention comparison. `pythia8_pdfdata_dir()` is the
compiled-in directory, overridden at run time by `$LIPOLGEN_PYTHIA8_PDFDATA`.

#### `DeuteronWaveSource` — AV18 from file, or CD-Bonn analytically

`DeuteronConvolutionB1::Options::wave` selects the deuteron wave function the
A = 2 gate is fed:

| value | what it is |
|---|---|
| `kFdeutFile` (**default**) | the tabulated `data/vmc/deuteron/fdeut.av18`, i.e. AV18, P_D = 5.76 % — **unchanged**, and a pytest asserts that setting this explicitly is a bit-for-bit no-op |
| `kCdBonn` | the analytic CD-Bonn parameterisation (Machleidt, PRC **63** (2001) 024001, Appendix D, Table XX) built by `cdbonn_wave()` / `cdbonn_fdeut_table()` in `cluster.hpp`; P_D = 4.86 %, Q_d = +0.2702 fm² against his published 4.85 % and 0.270 fm². `cdbonn_k_max_fm` / `cdbonn_dk_fm` (default 20.0 / 0.1) are `fdeut.av18`'s **own** grid, so an AV18 row and a CD-Bonn row differ in the wave function and in nothing else |

`C₁₁` and `D₉..D₁₁` are **computed** from the r → 0 boundary conditions, never
typed. CD-Bonn is CDKS's own wave function and closes the gate's remaining
factor 1.186 (above) — but it is **not** the default, because flipping it would
move every pinned gate number, and it is not a uniform improvement on an
arbitrary PDF.

### `--b1-unpol` — which unpolarised PDF the convolution folds against

`--b1-model li6-convolution` builds its own F₁ from an `UnpolSF`
(`Li6ConvolutionOptions::unpol`). Until 2026-09-04 that object was hard-wired
to the library `ToyF2` inside `default_inclusive_kernel`, so **every ⁶Li b₁
number the generator could emit was made on the toy** — including after the
A = 2 gate passed on CDKS's own MSTW2008 LO. `--b1-unpol` makes the other
backends selectable from the run surface, through the pipeline's own kernel
construction:

| value | backend | needs |
|---|---|---|
| `toy` (**default**) | the inclusive kernel's own `ToyF2`, handed to the convolution as **one shared object** — bit for bit what every published number was made with | — |
| `mstw` | `MstwSF()` — MSTW2008 LO over PYTHIA's `pdfdata` grid, the PDF CDKS computed their b₁ᵈ with | the optional PYTHIA tier |
| `ct18nlo` | `LhapdfSF("CT18NLO", 0)` — the phase-D stand-in, kept selectable so the PDF systematic can be quoted rather than remembered | the optional LHAPDF tier |

Measured on the shipped observable `Li6ConvolutionB1::b1(x, 2.5)` (2026-09-04):

| x | toy | mstw | mstw/toy | ct18nlo/toy |
|---|---|---|---|---|
| 0.10 | −2.796612e-05 | −5.167487e-05 | **1.847766** | 2.221302 |
| 0.30 | −1.391083e-04 | −1.774967e-04 | **1.275961** | 1.238833 |
| 0.50 | +1.034787e-04 | +8.453904e-05 | **0.816971** | 1.045870 |

Up to a factor 1.85 and **not monotone**, so the nucleon PDF is a shape change
in the ⁶Li observable and not a normalisation the 100 % band would absorb.

**Scope, and the price of it.** The flag reaches
`Li6ConvolutionOptions::unpol` and nothing else. `InclusiveKernel`'s own
`f2_source` — the F₁ of the unpolarised rate, and so the D_φ denominator of
the tensor weight — is set by a **separate** flag, `--unpol-sf` (§2b), and at
its shipped default `toy` it is `ToyF2` on every `--b1-unpol` setting, so the
**spin-blind rate is bit for bit** under *this* flag (`InclusiveSampler::cell_xsec_pb` is
array-identical, the same invariant the 100 % band satisfies) and what moves
is the tensor shift alone. Measured end to end on a ⁶Li inclusive
tensor-thirds run (seed 7, `--x-max 0.95`, default grid): the azz0 − azz±
cross-section split goes from **12.6037 pb** on `toy` to **10.6743 pb** on
`mstw`, a factor **0.8469**; `sigma_per_category_pb` and `sigma_pb` carry
that shift and are therefore not bit-identical. The consequence of the narrow
scope is that with anything but `toy` the numerator's F₁ and the denominator's
F₁ are no longer the same object, so the partial cancellation the `r_func`
default relies on is gone. Say which backend a plot used.

**And it is deliberately not merged into `--unpol-sf`.** The two are
different quantities with different conventions — this one is the
**deuteron's** per-nucleon F₁ inside CDKS Eq. (22) through `f1_cdks`, which is
explicitly *not* `NuclearF2::f1a`, while `--unpol-sf` is the whole-nucleus
F₂A/F₁A of the **rate** — and this one is a CDKS-comparability choice the
A = 2 gate verdict rests on (G3b = 0.843243 on MSTW, 0.440 on the toy). A
user who wants a realistic rate must not be forced to move the b₁ gate row off
MSTW as a side effect, nor the other way round. One refusal closes the
ambiguity the two flags create: `default_inclusive_kernel` folds the
convolution against the *kernel's* own `f2_source` when `b1_unpol` is null, so
under `--unpol-sf ct18nlo` a `--b1-unpol toy` would quietly mean CT18NLO while
`meta["b1_unpol"]` still wrote `"toy"`. `PipelineConfig::validate()` throws on
exactly that combination and asks for the b₁ backend to be named.

**It is never silently downgraded.** Selecting `mstw` in a build with no
PYTHIA tier, or with no `mstw2008lo.00.dat` on disk, raises at configuration
time and names what is missing — the tier, or the file and the directory it
was looked for in. The core library links neither optional tier, so the enum
(`PipelineConfig::b1_unpol`, the provenance `meta["b1_unpol"]` records) and
the object (`PipelineConfig::b1_unpol_sf`, the realisation) are two fields of
one choice, exactly as `optics_choice` / `optics` are; `set_b1_unpol` builds
the object, and `validate()` refuses a named backend with an empty slot.
Assigning `config.b1_unpol_sf` yourself names the choice `custom`.

`validate()` also refuses the flag with `--b1-model miller` and `cdks`, which
never read it — the same provenance rule as `--b1-band-scale` — and the
existing refusal of a **caller-supplied `kernel` together with a non-miller
`b1_model` is unchanged**: the point of this flag is that the gate-passing
configuration no longer needs a hand-built kernel to be expressed.

```bash
lipolgen-run --b1-model li6-convolution --b1-unpol mstw \
             --b1-band-scale 0 --x-max 0.95 --events 100000 --npz b0.npz
```

```python
cfg = lipolgen.make_config(b1_model="li6-convolution", b1_unpol="mstw",
                           events=100000)
```

### The top x cell — why `--x-max 0.95`

Both opt-in backends carry the CDKS camp's b₁ᵈ, which is a **Q² = 2.5
digitization with no Q² evolution**. In the topmost cell of the default window
(x = 0.954993) F₁ has fallen far enough that b₁/F₁ reaches **6.52** for `cdks`
and **5.91** for `li6-convolution` — past where the phi-averaged density 1 + w_avg
stays positive, and `InclusiveSampler` refuses the run with *"negative
phi-averaged density"*. Use `--x-max 0.95` (or `cfg.scenario.x_max = 0.95`).
`lipolgen-run` now checks this **before** building the pipeline and says which
flag fixes it, rather than letting that message — which names neither `--x-max`
nor the CDKS table — come out of the sampler three frames down.
`miller`'s b₁ is a ratio model (0.145 at the same point) and does not need it.
Those three figures are MEASURED off the shipped kernel (re-measured
2026-09-16) and are defined once, in `python/lipolgen/cli.py`'s
`B1_TOP_CELL_B1_OVER_F1`, which the `--x-max` help and the `--b1-model`
refusal both format their sentence from. Until 2026-09-16 this page said
6.6 / 5.6, that help said 6.6 / 5.6, the refusal said 3.3 / 5.6 and
`python/tests/test_b1_model.py` said 3.3 / 5.6 — the 3.3 being the
pre-2026-09-03 half of the `cdks` figure and the 6.6 twice its rounding, so
none of the four was what the code computes.
This is a property of the CDKS camp's table, not of the wiring: a hand-built
kernel with `Li6B1(CdksB1())` behaves identically.

### The four terms, and the numbers

The model (design §1.7), per nucleon:

> b₁^{⁶Li} = (2/6) ∫ (dz/z) { [f_S(z) + w_CG f_D(z)] b₁ᵈ(x/z) + w_αd δ_T f_αd(z) F₁ᵈ(x/z) }
> + (4/6) w_αd ∫ (dz/z) δ_T f_α(z) F₁^α(x/z)

with w_CG = **1/10 exactly** (the Clebsch–Gordan tensor dilution of the
deuteron inside the L = 2 α–d component) and w_αd = 1 nominal. The four terms
are reachable separately — `b1_embedded_s` (1), `b1_alpha_d_dwave_d` (2d),
`b1_alpha_d_dwave_alpha` (2α), `b1_cg_dwave` (3) — and they sum to `b1` to
1e-12. The α is J = 0 so b₁^α ≡ 0, but it is a **constituent** of the ⁶Li
spectral function whose light-cone density carries the same (3cos²θ − 1)
alignment, so term (2α) exists and is ≈ 0.5 × term (2d): a regression that
drops it moves b₁ by ~30 %.

The measured tables — the four terms at five x, the band rows, the knob rows,
`finite_q_delta` on/off, the densities, the ±5 % N_αd systematic and the T14
truncation remainders — are in
`docs/open_items/run_2026-09-02/phase_D_numbers.md` and are summarised in
`docs/OPEN_ITEMS_SOLUTIONS.md` §10.

**They WERE REGENERATED ON 2026-09-06 — as a BAND, which is the only form
they may be quoted in** (`docs/open_items/run_2026-09-06/phase_A_li6_tables.md`;
one row pinned at rtol 1e-12 by `python/tests/test_li6_unpol_band.py`). The
regeneration reran **the same observables at the same (x, Q²) points with the
same options** over `--b1-unpol {toy, ct18nlo, mstw}` × `--unpol-sf {toy,
ct18nlo, mstw}`, so the `toy` column of `phase_D_numbers.md` is reproduced to
every printed digit and the other two are a like-for-like band around it.

**Two bands ride with every ⁶Li b₁ number, and neither contains the other at
low x.** The unpolarised-backend envelope at Q² = 2.5 is

| x | envelope of x·b₁ over `--b1-unpol` | × {0, 1, 2} — what may be quoted | backend spread |
|---|---|---|---|
| 0.05 | [+1.574759e−6, +3.650716e−6] | **[0, +7.301431e−6]** | ×2.318 |
| 0.10 | [−6.212120e−6, −2.796612e−6] | **[−1.242424e−5, 0]** | ×2.221 |
| 0.20 | [−2.803165e−5, −2.120522e−5] | **[−5.606329e−5, 0]** | ×1.322 |
| 0.30 | [−5.324901e−5, −4.173248e−5] | **[−1.064980e−4, 0]** | ×1.276 |
| 0.50 | [+4.226952e−5, +5.411263e−5] | **[0, +1.082253e−4]** | ×1.280 |

and the **mandatory ±100 % A > 2 band** — which comes from
Q(⁶Li) = −0.0818(17) fm² against Q_d = +0.2859(3) fm², **not** from the A = 2
gate, and which the gate cannot lift because it is an A = 2 gate — is the
right-hand column. **At x = 0.05 and 0.10 the PDF spread exceeds ×2, so
running `{0, 1, 2}` on one backend does not bracket the others there**; at
x = 0.20/0.30/0.50 it does. Quote **both**, and never one row.

**Which axis moves what, measured.** `--b1-unpol` moves terms (2d) and (2α)
and nothing else — terms (1) and (3), the densities, the T14 remainders and
the `--b1-alpha-d-dwave-weight 0` row are **bit-identical** across it, because
they convolve b₁ᵈ rather than F₁. `--unpol-sf` is **exactly flat on x·b₁**
(`Li6ConvolutionB1::b1` ignores its `f1` argument) and moves the *denominator*
of A_zz instead: A_zz itself spans a factor **2.23** over the seven legal
cells at x = 0.10, Q² = 2.5, y = 0.5. Two of the nine cells are refused by
name (`--b1-unpol toy` under a non-toy `--unpol-sf`; the row is in §2b's
"What is refused, and why"), so the band on these tables is **three wide, not
seven**.

**Two published readings did NOT survive the regeneration and are withdrawn as
readings** (both are flagged in place in `phase_D_numbers.md` and
`OPEN_ITEMS_SOLUTIONS.md` §10): *"terms (2d)+(2α) have the opposite sign to
term (1) over the whole window"* is a `toy` feature — on `ct18nlo` the orbital
sector has the **same** sign as term (1) at every x, and on `mstw` at four of
five — and *"w = 2 nearly zeroes b₁ at x = 0.10"* is toy-only, since on the
two real fits `--b1-alpha-d-dwave-weight 2` moves b₁(0.10) **away** from zero
by +24 % and +9 %. `∫b₁ dx` over [0.01, 1.2] becomes
**[+1.387049e−4, +1.458804e−4]** (a 5.17 % spread, the tightest quantity on
the axis).

**Why the rerun was a rerun and not a factor**, which is the part that has not
changed: they were recorded under a publication ban that the A = 2 gate lifted
on 2026-09-03 — but the lift is a statement about a *configuration*, and the
`toy` column is still not in it. Three facts decide that, and they are why
"quotable with `r_sigma_lt` and `finite_q_delta = false`" was not enough:
naming the two small knobs and omitting the unpolarised backend named the
wrong things.

1. Every row was made with the **default `ToyF2`** unpolarised input. That
   configuration's own G3b is **0.440** — *outside* the [0.5, 2] acceptance
   window the lift was read off. Nothing that fails the gate's own magnitude
   clause is covered by a lift granted on a row that passes it.
2. The gap is a **shape change, not a normalisation**: mstw/toy on
   `Li6ConvolutionB1::b1(x, 2.5)` is **1.848 / 1.276 / 0.817** at
   x = 0.10 / 0.30 / 0.50. No single factor converts the ToyF2 tables into
   `mstw` ones, so they cannot be rescaled and must be recomputed. Checklist
   item 4 puts the same effect at **×1.92** on the gate's own peak.
3. What survives regeneration untouched, and may be cited as it stands, is
   only what is **algebraic** and so configuration-independent: the four terms
   summing to `b1` to 1e−12, and the band's exact linearity (the scale-2 row
   is bit-for-bit 2× the scale-1 row). Every **value**, every ratio *between*
   terms, and every derived systematic in those tables is a ToyF2 measurement
   and is labelled as such until it is rerun.

The standing rules apply to the regenerated tables exactly as they did before:
**always as a {0, 1, 2} × b₁ band and always with the configuration that made
them** — which means naming the unpolarised backend first, and then
`r_sigma_lt` and `finite_q_delta = false`, the ⁶Li backend's defaults and not
the gate's. What the regeneration adds is that the backend is now a **band and
not a single name**, and that the two bands must both be quoted.

## 2b. Structure-function backends for **every** channel — `--unpol-sf` / `--pol-sf`

```
lipolgen-run --channel inclusive        --unpol-sf ct18nlo --pol-sf nnpdfpol
lipolgen-run --channel tagged-6Li-alpha --unpol-sf ct18nlo
```

```python
cfg = lipolgen.make_config(channel="tagged-6Li-alpha", events=100_000,
                           unpol_sf="ct18nlo", pol_sf="nnpdfpol")
```

**What they are.** `--unpol-sf {toy,mstw,ct18nlo}` chooses the unpolarised
backend that supplies F₂ — and through it F₁, F_L and the whole unpolarised
rate — and `--pol-sf {toy,nnpdfpol}` the polarised backend that supplies g₁,
and through the Wandzura–Wilczek relation g₂. Both default to the toy backends
and are then **bit for bit** what every published number was made with.

**They do not have the same reach, and the run says so.** `--unpol-sf` reaches
**every kernel the pipeline builds**: the inclusive kernel, the coherent
channel (which owns no structure function of its own and rides the inclusive
cell cross sections — ×0.705841 on the shipped ⁶Li `ct18nlo` coherent run),
and the tagged struck-cluster kernel — `read` on every channel and every plan
in the knob-provenance table. `--pol-sf` reaches the inclusive and the tagged
kernels **only where the fill also carries `lam_e · P_e ≠ 0`**, and not the
coherent channel at all, whose rate is spin-independent. Both halves matter:
under this CLI's own default plan (`tensor-thirds`, every category at
`lam_e = 0`) `--pol-sf` is read on NO channel. See "`--pol-sf` does not run on
the coherent channel, nor under an unpolarised-beam plan" below.

**Why they exist.** Until 2026-09-04 the run surface could not select either.
`default_inclusive_kernel` hard-wired `ToyF2`, and `struck_cluster_kernel` had
*no* structure-function slot of any kind — so a tagged run was `ToyF2`/`ToyG1`
with no way to say otherwise, and not even the hand-built
`PipelineConfig::kernel` escape hatch reached it (the `Pipeline` reads that
field on the non-tagged branch alone; §7 below). `sf.hpp` labels both
defaults **TOY** and anchors `ToyF2` *by eye* — "adequate for phase-space maps
and factor-1.5 rate estimates ONLY".

### What the toy has been costing, measured

⁶Li, `default_configs("6Li")[1]` (e 10 GeV × ⁶Li 99.5 GeV/u), the shipped
window, run through the pipeline on 2026-09-04. σ is the accepted cross
section `sigma_pb()`; A_zz is the canonical thirds estimator on
`tensor_thirds_plan(0, 0.6)`'s per-category cross sections divided by P_zz;
A_∥ is `(σ₊ − σ₋)/(σ₊ + σ₋)/(P_e P_z)` on `helicity_flip_plan(1, 0.7, 0.7)`.
**Every column is at the shipped `Scenario::x_max` = 1.0** — the A_zz column
was printed at `--x-max 0.95` until 2026-09-05, in rows whose σ and A_∥ were at
1.0, so one row carried two windows (at 0.95 the three A_zz are
−5.193192e−4 / −6.503956e−4 / −6.545724e−4, and both ratios are unchanged to
six digits, which is what hid it).

| channel | backend | σ [pb] | ratio | A_zz | ratio | A_∥ | ratio |
|---|---|---|---|---|---|---|---|
| inclusive | `toy` | 591846.2 | 1 | −5.193231e−4 | 1 | −1.171716e−3 | 1 |
| inclusive | `--unpol-sf ct18nlo` | 472571.9 | **0.79847** | −6.503970e−4 | **1.25239** | −9.668880e−4 | 0.82519 |
| inclusive | `--unpol-sf mstw` | 469556.5 | **0.79338** | −6.545739e−4 | **1.26044** | −1.057030e−3 | 0.90212 |
| inclusive | `--pol-sf nnpdfpol` | 591846.2 | 1 (exactly) | −5.193231e−4 | 1 (exactly) | −1.270806e−4 | **0.10846** |
| inclusive | both | 472571.9 | 0.79847 | −6.503970e−4 | 1.25239 | −1.591593e−4 | 0.13583 |
| tagged-6Li-α | `toy` | 591846.2 | 1 | **0 exactly** | — | −3.514730e−3 | 1 |
| tagged-6Li-α | `--unpol-sf ct18nlo` | 472571.9 | **0.79847** | **0 exactly** | — | −2.900243e−3 | 0.82517 |
| tagged-6Li-α | `--unpol-sf mstw` | 469556.5 | **0.79338** | **−1.3774e−16** (qual. 1) | — | −3.170626e−3 | 0.90210 |
| tagged-6Li-α | `--pol-sf nnpdfpol` | 591846.2 | 1 (exactly) | **0 exactly** | — | −3.811967e−4 | **0.10846** |

Read the three qualifications with the table.

1. **The tagged A_zz is zero on every backend by construction, not by
   accident — and on `mstw` it is a floating-point zero, not a bit-level
   one.** At the shipped default the struck-cluster kernel has no b₁
   (`--inclusive-b1` is off, because on a tagged channel the α–d density is
   already in the event weight), so the *total* per-category cross sections
   carry no tensor term at all. On `toy`, `ct18nlo` and `nnpdfpol` the three
   totals are then **equal to the last bit** and the thirds estimator returns
   **0.0 exactly**. On `mstw` they are not: measured 2026-09-05 on
   `tensor_thirds_plan(0, 0.6)`'s `sigma_per_category_pb()`, they are
   469556.45134814945 / 469556.45134814945 / 469556.4513481495 — σ₀ is **one
   ULP** (5.82e−11 pb) above σ_±, so A_zz comes back **−1.3774e−16** rather
   than 0.0. It is round-off in the summation order — 2.7e−13 of the
   inclusive channel's own toy A_zz (−5.193231e−4) — and *not* a tensor
   signal. Quote "0 by construction", never "0 to the last bit",
   unless the backend is one of the three where the last bit is what was
   measured. The tagged tensor signal lives in the
   spectator-momentum-differential rate, not in σ_tot. Turn the inclusive b₁
   on (`--inclusive-b1 --x-max 0.95`) and it moves by the same factors the
   inclusive channel does: A_zz = −1.557969e−3 (`toy`) → −1.951190e−3
   (`ct18nlo`, ×1.25239) → −1.963721e−3 (`mstw`, ×1.26044).
   **`--inclusive-b1` and a TILTED plan are refused together on
   `tagged-6Li-alpha`** (since 2026-09-16): the struck cluster's |S_c m_S⟩ is
   evaluated with its quantization axis along the **beam** whatever the ion
   fill's axis is, so `--plan transverse-tensor|tensor-flip --inclusive-b1`
   used to run — and to record `inclusive_b1` as READ — while applying the
   struck deuteron's b₁ with P₂(cos 0) = +1 on a fill sitting at
   P₂(cos 90°) = −½. The refusal names the three ways out: run the tilted fill
   spin-blind, run the b₁ on an untilted fill, or use the inclusive channel,
   whose kernel does carry the axis. The same refusal catches P_e ≠ 0 on a
   tilted fill, on **all three** tagged channels. It does **not** fire for
   `--inclusive-b1` on the ⁷Li alpha tag (C++ `tagged-7Li-alpha`) or on
   `tagged-d-p`, where the struck cluster is spin ½ and the knob is measurably
   inert. **Name the RESOLVED channel here, never the CLI alias**:
   `--channel tagged-alpha` resolves against `--isotope` (the `CHANNELS`
   table in `python/lipolgen/__init__.py`), so on the default ⁶Li it IS
   `tagged-6Li-alpha` and `lipolgen-run --channel tagged-alpha --plan
   transverse-tensor|tensor-flip --inclusive-b1` exits **1** with this refusal
   (measured 2026-09-16, both tilted plans, `--events 2`). The same command
   line under `--isotope 7Li` also exits 1 — but on a different refusal,
   reached before this one: both shipped tilted plans are spin-1 patterns and
   ⁷Li is J = 3/2, so there is no tilted ⁷Li fill for this test to see. With
   `--channel tagged-d-p` it exits **0** under either isotope spelling (the
   channel implies the deuteron). Nothing spin-blind and
   nothing untilted moved: every CLI default, every reference gate and both
   shipped tilted plans (which carry P_e = 0) are bit-identical.
   **Why `--x-max 0.95` is in that command line** (measured 2026-09-05, and it
   was prescribed here without a reason until then): with `--inclusive-b1` on
   a tagged channel *and* a non-toy unpolarised backend, the shipped window's
   top cell (x = 0.955, Q² = 167.3) gives **1 + w_avg = −0.1302** on
   `ct18nlo` and **−0.03496** on `mstw`, and `InclusiveSampler`'s
   configuration-time positivity check refuses the run by name — the same
   refusal `--b1-model cdks|li6-convolution` already carries on the inclusive
   channel, reached here through a different door. On `--unpol-sf toy` the
   same command runs at x_max = 1.0, which is why the flag carries no
   unconditional refusal; the numbers above are all at 0.95 so that one
   window covers the four backends.
2. **The A_∥ ratios are window integrals of a SIGN-CHANGING integrand, and
   are smaller than any of their own parts.** On the toy the cross-section-
   weighted ⟨A_∥⟩ is −2.34618e−3 for x < 0.01 (53.0 % of the rate) and
   *positive* above x = 0.05, so the window number is a near-cancellation.
   Per x band the polarised swap is ×0.163 (x < 0.01), ×1.041
   ([0.01, 0.05)), ×0.905 ([0.05, 0.1)), ×1.172 ([0.1, 0.3)), ×0.937
   (x ≥ 0.3) — the window-integrated **0.108 is smaller than every band
   ratio**, because the bands partly cancel. Quote the band, not the ×0.108,
   unless the window is the observable.
3. **The two selectors are NOT orthogonal.** `InclusiveKernel` builds its
   default `ToyG1` on its *own* base `UnpolSF`, so `--unpol-sf` alone moves
   g₁ as well: at `--pol-sf toy`, `--unpol-sf ct18nlo` moves ⁶Li F₁ by
   ×0.9942 / ×1.0711 / ×1.2038 / ×1.0014 and g₁ by ×1.0544 / ×1.1346 /
   ×1.3050 / ×1.0551 at x = 0.05 / 0.10 / 0.30 / 0.50, so A₁ = g₁/F₁ moves
   5–8 %. **`--pol-sf toy` does not mean "g₁ unchanged".** State both
   settings next to any number.

### `--pol-sf` does not run on the coherent channel, nor under an unpolarised-beam plan

**Two axes, and the second one is the shipped default.** g₁ enters the rate
through exactly one product — `InclusiveKernel::amplitudes` adds
`lam_e · P_e · (m/J) · cos θ_S · A_∥` and nothing else in a run reads a
`PolSF` — so the selector is read only where the CHANNEL evaluates g₁ **and**
the FILL carries `lam_e · P_e ≠ 0`.

*The channel.* The coherent yield is `f_coh(x)` times the **unpolarised** cell
cross sections, and that channel's tensor signal is the recoil azimuth's
1 + c₂ cos 2(φ_t − φ_S), which `CoherentSampler` owns; not even `RcModel` asks,
its `applies()` being false there. Measured 2026-09-05 (`--channel coherent
--events 400 --seed 11`): `sigma_pb`, `sigma_per_category_pb` and **all 47
generated array columns** are bit-identical between `--pol-sf toy` and
`--pol-sf nnpdfpol`.

*The run plan.* `tensor_thirds_plan`, `transverse_tensor_plan` and
`tensor_flip_plan` build **every** category at `lam_e = 0, pe = 0`
(`src/core/bookkeeping.cpp`), and `helicity_flip_plan` at `--pe 0` is the
same — so on those plans g₁ is multiplied by zero on every channel.
`tensor-thirds` is this CLI's **default plan**. Measured 2026-09-05, 600
events seed 7, sha256 over all 47 columns plus both σ vectors: `--pol-sf
nnpdfpol` is bit-identical to `toy` on inclusive-⁶Li, inclusive-d,
tagged-⁶Li-α and tagged-d-p under `tensor-thirds`, and on inclusive-⁶Li under
`transverse-tensor`, `tensor-flip` and `helicity-flip --pe 0`. It **does**
move under `helicity-flip` at `--pe 0.7` on every channel but the coherent
one.

Neither axis is **refused**: `--unpol-sf` reaches the coherent rate, so
refusing its partner on one channel of a three-channel scan costs more than it
buys — and refusing it on the plan axis would make `--pol-sf nnpdfpol` fail at
the CLI's own defaults. Both are **labelled**, and §7c is where the label
lives now:

* `meta["pol_sf"]` carries `not read on channel coherent-6Li` or `not read
  under this run's fill` instead of a backend name — and that is what a
  `--pol-sf toy` run records too, because `ToyG1` did not run there either;
* `meta["pol_sf_reach"]` carries the whole sentence, from the one definition
  (`pol_sf_reach_report`);
* the banner's `KNOB PROVENANCE` block prints the same row, with the reason
  spelled out whenever the value is off its default.

The kernel still *carries* the selected `g1_model`, so
`p.dis_sampler.kernel.tables(x, q2).g1` does move; what never happens is that
the run reads it. The predicate is `pol_sf_is_read(config, plan)` — the plan
is an argument, and that is round 2 of §D6.

**`--pe` is a different question, with a different answer.** Under the three
tensor plans the flag is not read at all (their factories never take it), so
its row is labelled too; under `helicity-flip --pe 0` the plan *did* consult
it — that is what made the product zero — so `pe` is `read` there while
`pol_sf` is not.

### The neutron is a sign, not a factor

`ToyG1`'s a1n(x) = −0.07(1−x)² + 0.8x^2.2 crosses zero near x ≈ 0.25 and is
**positive** above it, while NNPDFpol1.1's g₁ⁿ stays negative to x ≈ 0.6.
Measured at Q² = 10:

| x | `toy` g₁ⁿ | `nnpdfpol` g₁ⁿ | ratio |
|---|---|---|---|
| 0.05 | −0.242789 | −0.248012 | 1.0215 |
| 0.10 | −0.0800273 | −0.136452 | 1.7051 |
| 0.20 | −0.0113549 | −0.0607773 | 5.3525 |
| 0.30 | **+0.00520678** | **−0.0273959** | −5.2616 |
| 0.50 | +0.00778817 | −0.000370697 | −0.0476 |
| 0.70 | +0.00247092 | +0.00221761 | 0.8975 |

So **the shipped toy g₁ⁿ has the wrong sign over roughly 0.25 < x < 0.6.** On
isoscalar ⁶Li the proton term dominates the sum and this mostly hides; on a
**neutron-tagged** run (`--channel tagged-d-p`, whose `dis_target` is
`NEUTRON_TARGET`) it does not. That is the single strongest reason `--pol-sf`
exists.

#### …and on that channel the pair needs `--x-max 0.95`

```
lipolgen-run --isotope d --channel tagged-d-p --plan helicity-flip \
             --pz 0.7 --pe 0.7 --unpol-sf ct18nlo --pol-sf nnpdfpol \
             --x-max 0.95
```

Without `--x-max 0.95` that command — the A_∥ measurement on a free neutron,
at the CLI's own default `--pe 0.7` — is **refused at configuration time**,
before an event is drawn. It is not a bug in either backend; it is what the
two of them do together at the top of the shipped window. Measured on
2026-09-05 in the cell the sampler names, x = 0.954993, Q² = 1119.1, on the
`tagged-d-p` DIS target (a free neutron):

| backend pair | F₁ⁿ | g₁ⁿ | A₁ = g₁/F₁ | 1 + w_avg at P_z = P_e = 0.7 |
|---|---|---|---|---|
| `toy` + `toy` | 1.37602e−5 | 9.94559e−6 | 0.722778 | builds |
| `ct18nlo` + `toy` | 1.19311e−6 | 8.62353e−7 | 0.722778 | builds |
| `toy` + `nnpdfpol` | 1.37602e−5 | 5.25659e−6 | 0.382014 | builds |
| **`ct18nlo` + `nnpdfpol`** | 1.19311e−6 | 5.25659e−6 | **4.40579** | **−0.02414** |
| **`mstw` + `nnpdfpol`** | 1.10033e−6 | 5.25659e−6 | **4.77729** | **−0.1105** |

`ToyG1` is A₁(x)·F₁ on the kernel's *own* F₁, so its A₁ is the same 0.722778
whatever supplies F₂ — the ratio only runs away when a **grid** g₁ is divided
by a **different** grid's F₁, which is exactly what the pair does. The
phi-averaged density 1 + w_avg is then negative and `InclusiveSampler` refuses
the run rather than draw max(W, 0), which would dilute the modulation *and*
skew the (x, Q²) mixture.

There is **no clamp**, and there will not be one: clamping g₁ or the density
is a physics change, and a silent one. What the refusal does instead is name
the cell, the A₁ that made it negative, and the cure —

```
negative phi-averaged density for m=-1 at x = 0.955, Q2 = 1119
(1 + w_avg = -0.02414).  That cell's A1 = g1/F1 is 4.406 (F1 = 1.193e-06,
g1 = 5.257e-06) at lam_e = 1, P_e = 0.7, J = 1: ...  Cure: lower the
acceptance window's x_max below 0.955 (Scenario::x_max; the CLI's --x-max,
e.g. --x-max 0.95 on the shipped grid, whose top cell is x = 0.955), or
reduce P_e / P_z, which scale w_avg linearly.
```

— and `lipolgen-run` prints that message and exits 1 as a **configuration
error, not a traceback** (the `PythiaBridge` rule: show the message the model
wrote). `--pe 0` also removes it, because w_avg scales with P_e·P_z; so does
either backend alone. It is the same class of refusal `--b1-model
cdks|li6-convolution` already carries in the *tensor* sector (§2a, "The top x
cell"), in the vector sector and on a different channel.

### The grid clause — `--unpol-sf ct18nlo` runs a third of the window off-grid

CT18NLO's grid starts at Q = 1.295 GeV, i.e. **Q² = 1.677**, and the shipped
⁶Li window's accepted cells start at Q² = 1.054. Below its grid LHAPDF does
**not** freeze — it continues the evolution downward and F₂ᵖ falls fast (at
x = 3e−4 it is a factor **2.34** below the toy at Q² = 0.7, and 1.70 below at
Q² = 1.0). Measured on the shipped window: **36.18 %** of a CT18NLO
**inclusive** run's own accepted rate (42.33 % of the *toy* run's) sits below
that floor.

**The fraction is per channel.** Its denominator is
`Pipeline::cell_rate_weights_pb()` — *this run's* per-cell rate — and on the
coherent channel that is σ_cell·f_coh(x) over the cells that admit a
diffractive mass, not the inclusive cells. `f_coh` falls by a factor ≈26
across the window, so the reweighting is large: measured 2026-09-05 on
`--channel coherent --unpol-sf ct18nlo` at config 1, **44.746 %** of the
coherent rate is below the floor against 36.179 % of the inclusive cells it
rides on, and those coherent weights sum to 8806.096207 pb, which is that
run's `sigma_pb()` exactly. Until 2026-09-05 a coherent run recorded 36.18 %,
a quantified `meta` key that did not describe the run it was attached to.

It is **not refused** — refusing would make the flag unusable on the shipped
scenario — but it is never silent either:

* the run banner prints the floor and the fraction on every run;
* `meta["unpol_sf_grid_q2_min"]` and `meta["unpol_sf_below_grid_frac"]` carry
  them into the npz;
* the fraction is computed from the loaded grid's own `q2Min()`
  (`LhapdfSF::q2_min`), never from a typed-in number, and against **this
  channel's** own per-cell rate (`Pipeline::cell_rate_weights_pb`), never
  another channel's;
* a backend that reports **no** floor writes **NaN**, not 0 — PYTHIA's
  `MSTWpdf` keeps its `qsqmin` private, and a `Custom` object reports nothing.

The two selectors' other `meta` keys are `unpol_sf`, `pol_sf` and
`pol_sf_reach` — the last two carrying the label and the reason on a channel
where g₁ was never evaluated, exactly as `b1_model` / `rank2_input` do for a
rank-2 slot nothing filled.

**There is no `--q2-min` flag** — `--x-max` is the only scenario knob on the
command line — so moving the window off the extrapolation means raising
`Scenario::q2_min` from Python (`sc = cfg.scenario; sc.q2_min = 1.7;
cfg.scenario = sc`), and it costs rate: on the shipped ⁶Li `ct18nlo` run
q2_min 0.7 → 1.7 takes the below-grid fraction to exactly **0** and σ from
472571.9 pb to **301599.8 pb** (×0.638). MSTW's own grid
boundary has **not** been characterised the way CT18NLO's has; its below-grid
F₂ᵖ was probed at two points (0.492 and 0.550 at x = 3e−4, Q² = 0.7 and 1.0
against the toy's 0.681) and looked milder, which is an observation and not a
measured fraction.

### What is refused, and why

| refusal | reason |
|---|---|
| a named backend with an empty object slot (`mstw`, `ct18nlo`, `nnpdfpol`, `Custom`) | the core library links neither the PYTHIA nor the LHAPDF tier and cannot build the backend; falling back to the toy would put the toy's numbers under a label that says otherwise. The message **names the missing tier**, and the CLI refuses one layer earlier still (`require_unpol_sf_tier` / `require_pol_sf_tier`) |
| `toy` with an object attached | the provenance and the realisation drifted — the attached backend would be silently dropped. Assigning `config.unpol_sf_obj` from Python sets `Custom` for you |
| either selector together with a caller-supplied `PipelineConfig::kernel` | the kernel wins on the inclusive branch and is not read at all on a tagged one; `meta` would record a backend that did not run |
| a directly set `struck.f2_source` / `struck.g1_model` with the matching selector still at `toy` | a C++ caller may put an object straight into the tagged kernel's slots (`Pipeline` fills them only when empty, the `breakup.triton_sf` arrangement), but `meta` would then record `"toy"` for a run made on something else. Name the selector `Custom` alongside |
| a caller-supplied `PipelineConfig::kernel` on a **tagged** channel, with or without a selector | the tagged channels draw from the **struck-cluster** sampler and never read that field, while `meta["unpol_sf"]`, `["pol_sf"]` and `["b1_model"]` all wrote `"caller-supplied kernel"`. Measured: a `tagged-6Li-alpha` config carrying an `InclusiveKernel` built on CT18NLO ran at F2A(0.3, 10) = 0.3691149345, the **toy** value. Still accepted on `inclusive` **and `coherent`**, which do read it |
| a `PythiaBridge` whose `options.f2_source` is not the config's own `unpol_sf_obj` (`set_pythia_hadronizer`) | the T2 struck-nucleon species draw would be on a backend `meta["unpol_sf"]` does not name — up to 24 % on F₂ⁿ/F₂ᵖ. Identity, not equality: hand the bridge `config.unpol_sf_obj` before constructing it, which `lipolgen-run` and `lipolgen.run` do |
| `--b1-unpol toy` under a non-toy `--unpol-sf`, on `--b1-model li6-convolution` | `toy` there means "the kernel's own `UnpolSF`, shared as one object", and that object is no longer `ToyF2` — so `meta["b1_unpol"]` would record `"toy"` for a b₁ folded against CT18NLO. Name the b₁ backend explicitly (§2a) |

What is **not** refused: either selector on a **tagged** or the **coherent**
channel. Unlike `--b1-model` — which `validate()` refuses off the inclusive
channel and off ⁶Li — `--unpol-sf` reaches the rate on every channel; that is
the whole point of it, and the provenance table says `read` for it on every
channel and plan. `--pol-sf` does **not** reach the coherent rate, and does
not reach any channel under an unpolarised-beam plan (`tensor-thirds`, this
CLI's default, included). It is not refused on either axis: it is *labelled*,
so `meta` records which of the two reasons applies instead of naming a
backend. A rule that depends on the run PLAN could not be refused in any case
— `validate()` has no plan. See "`--pol-sf` does not run on the coherent
channel, nor under an unpolarised-beam plan" above.

Also not refused: the neutron-tagged `ct18nlo`/`mstw` + `nnpdfpol` pair. That
one is refused by the **sampler**, not by `validate()`, and only where it
actually fails — the acceptance window's top cell — with the cell and
`--x-max` in the message. See "…and on that channel the pair needs
`--x-max 0.95`" above.

### Scope, and the three things this does not cover

* **R stays its own axis.** Nothing in the pipeline sets an `RFunc`, so R is
  `r_sigma_lt` at all four places `InclusiveKernel` needs it, on every
  channel and every backend. It is not folded into these two flags because
  it is not a small correction. Measured over the 3051 accepted cells of the
  shipped ⁶Li window on 2026-09-04: `r_sigma_lt` spans 0.0046–0.1763 and
  `r1998` spans 0.0124–0.4010, the ratio (1 + r1998)/(1 + r_sigma_lt) spans
  **0.9203–1.1910** with a cross-section-weighted mean of **1.1229**, and
  **38.18 %** of the accepted cell cross section lies outside R1998's own
  stated support (`R1998_X_MIN/MAX` = 0.005/0.86, `R1998_Q2_MIN/MAX` =
  0.5/130), where `r1998(..., clip = true)` returns the clipped boundary
  value. So R1998 is not a drop-in default, and shipping an `--r-model` in
  the same change as these two would make three independent 10–40 %
  movements unattributable. §2a, "Which R each object uses", has the
  decision that keeps `Li6ConvolutionOptions::r_func` on `r_sigma_lt`.
* **The EMC hook is still empty in every kernel.** `Options::emc_ratio` is
  set by no site in `src/`, `python/` or `examples/`, so no shipped run
  applies any medium modification to F₂A on any channel. If one is ever
  wired, note that the transcribed EPPS21 depletion constant
  (`EMC_VALENCE_DEPLETION_EPPS21`) is quoted against **CT18ANLO**, which is
  *not installed here* — while `Epps21Ratio`'s own `proton_set` default is
  CT18NLO. Referencing a ratio and a baseline to two different proton fits is
  the same 4.2 %-class error that constant's provenance note already
  describes, so any EMC hook must state its own proton denominator.
* **`--b1-unpol` is a separate flag** and stays one; §2a says why, and the
  one refusal that keeps the two from mislabelling each other.

### One backend for the whole run, including T2

The T1 breakup's struck-nucleon species draw takes the sampler kernel's own
`UnpolSF` object (`BreakupOptions::f2`, filled by the `Pipeline`), and
`InclusiveGenerator::proton_fraction` and `InclusiveSampler`'s cell cross
sections read the same object — so all three follow the selector for free and
must not grow knobs of their own. The **T2 PYTHIA bridge** is the one consumer
outside the `Pipeline`'s reach: `lipolgen-run` builds the bridge itself, so it
now hands `PythiaBridgeOptions::f2_source` the same object it gave the kernel.
Without that a `--hadronize --unpol-sf ct18nlo` run would draw its T2 species
from `ToyF2` while its rate came from CT18NLO, and the toy's F₂ⁿ/F₂ᵖ — the
straight line clip(1 − 0.75x, 0.25, 1) — is 0.9625 / 0.8500 / 0.6250 at
x = 0.05 / 0.20 / 0.50 against CT18NLO's 0.9218 / 0.7219 / 0.5035, i.e. up to
**24 %** away. `lipolgen.run(hadronize=True, ...)` does the same wiring (an
explicit `pythia_options={"f2_source": ...}` still wins *there*).

A caller building a bridge by hand no longer *may* forget:
`set_pythia_hadronizer(config, bridge)` **refuses** a bridge whose
`options.f2_source` is not the config's own `unpol_sf_obj`, naming both sides
and the fix. It is object **identity**, the rule everything else in a run
follows (the `li6-convolution` b₁ shares the kernel's own `UnpolSF`;
`BreakupOptions::f2` is the sampler kernel's own) — two separately constructed
`LhapdfSF("CT18NLO", 0)` are two grids whose agreement nothing checks. Both
unset is the same object, so the default run is untouched.


## 2c. ⁷Li inclusive — the tensor sector of that run is **exactly zero**

Read this before you plan a ⁷Li A_zz or cos 2φ measurement.

```
$ lipolgen-run --isotope 7Li --plan helicity-flip --events 100
  7Li rank-2: EMPTY -- 7Li is spin 3/2 and default_inclusive_kernel fills a
     rank-2 slot for spin 1 ONLY, so b1_32 = b2_32 = delta_32 = 0 and the
     tensor term of the phi-averaged rate, the cos 2phi (gluon transversity)
     amplitude and therefore A_zz of this run are IDENTICALLY ZERO, not small
     ...
```

**What is zero, and why.** `InclusiveKernel::tables` dispatches the rank-2
slots on the ion spin — spin 1 reads `b1_func / b2_func / delta_func`, spin 3/2
reads `b1_32_func / b2_32_func / delta_32_func` — and an unset slot is `0.0`
(an unset `b2` is 2x·b₁ and therefore 0 too). `default_inclusive_kernel` fills
the **spin-1** slots only, because **there is no published b₁ for ⁷Li**, so on
a ⁷Li inclusive run the tensor term of the φ-averaged rate, the cos 2φ
amplitude and A_zz are all identically zero. Measured on the honest ⁷Li A_zz
plan (pure |m| = 3/2, T = +1, against pure |m| = ½, T = −1; 60 000 events,
seed 11, config 1): the two per-category cross sections come back as **the same
double**, 590952.42641509 pb, and the asymmetry is **exactly 0.0**. There is no
"small" here — the structure functions are absent, not suppressed.

**And that is why the ⁷Li fill's alignment reaches nothing.** Measured
2026-09-06 (`docs/open_items/run_2026-09-06/phase_A_numbers.md` §A3): moving
this run's fill from the max-entropy ladder (T = 0.4) to an explicit
`--pzz-mode typed --pzz 0.5` changes **no observable** — σ agrees to **1 ulp**
(1.970 × 10⁻¹⁶ relative), one per-category σ exactly and the other to
3.9 × 10⁻¹⁶, and A_∥ moves by 6.2 × 10⁻¹⁴ of its own statistical error —
because the only fill moments the kernel sums are Σ p_m = 1 and ⟨J_z⟩/J = P_z,
and **both fills honour `--pz`**. The **sample** is not identical, though:
0.1060 % of 100 000 events (106) land in a different (x, Q²) cell, since cell
selection is a discrete function of weights that differ in the last bit. What
the typed fill does buy is the **recorded T**, 0.4 → 0.5, the divisor of any
tensor estimator — worth −20 % on δ(A_zz) and δ(cos 2φ) at fixed N, on an
asymmetry that is identically zero here.

**Nothing else about the run is affected.** The unpolarised rate and the whole
vector sector (g₁, A_∥) are correct at J = 3/2, and so is every piece of
machinery around the missing input: with a rank-2 slot supplied by hand on a
caller-built kernel (`b1_32_func = +0.05·F1`, `delta_32_func = −1e−2·F1`,
40 000 events, seed 11) the sampler returns σ = 565389 / 616516 pb and
**A_T = −0.043258**, which is −b₁/F₁ = −0.05 times the sample's own
1/(1 + εR) ≈ 0.865 (`docs/theory/SPIN32_FINITE_GAMMA.md` Eq. (48)). What is
missing is a physics **input**, not code.

**How the run says so** (it used to say nothing at all, while `meta` recorded
`b1_model = "miller"` — a backend that did not run):

* the banner block above, on **every** ⁷Li inclusive run, tensor plan or not;
* `meta["rank2_input"]` — the same sentence, from the same C++ definition
  (`rank2_input_report`), on every run and every channel;
* `meta["b1_model"] = meta["b1_unpol"] = "none (spin 3/2: no rank-2 input)"`.
  The two **scales** stay numeric at 1.0 — `validate()` already refuses any
  other value on that path, so they cannot record a variation that did not run.

**What to do instead.**

* **An inclusive tensor programme needs `--isotope 6Li`**, where `--b1-model`
  reaches a real backend (§2a) — with its band, and with the configuration
  quoted.
* **A ⁷Li tagged run is not affected.** `--channel tagged-7Li-alpha` carries
  the α–t alignment in the event weight, and it is gated today at
  ⟨P₂(cos θ_k)⟩ = −T/5. `meta["rank2_input"]` says so in its own words there
  rather than borrowing this sentence; so does the coherent channel, whose
  tensor signal is the recoil azimuth.
* **There is no spin-3/2 tensor plan in this tree.** `--plan tensor-thirds`,
  `transverse-tensor` and `tensor-flip` are spin-1 patterns (their categories
  hard-code j = 1) and are **refused** at J = 3/2 — two of them used to build a
  spin-1 plan that died inside the sampler two frames later. `helicity-flip` is
  the only plan that takes `j`, and it is a **vector** plan. The honest ⁷Li
  A_zz plan is the two-state T = +1 / T = −1 contrast, built by hand.
* **`--pzz` is read by the three spin-1 tensor plans always, and by
  `helicity-flip` only at `--pzz-mode typed`.** At the default
  `--pzz-mode ladder` that plan builds its fill from the max-entropy ladder at
  `--pz`, so `--plan helicity-flip --pzz 0.6` gives T = 0.4 at J = 3/2, not
  0.6. The banner prints the fill's own moments — and, on `helicity-flip`, the
  fill the *other* mode would have built — on every run, so the difference is
  never silent. `--pzz-mode typed` (2026-09-06) honours the typed value and
  **refuses**, with the edge named, one outside the plan's domain rather than
  clamping to it: at `--pz 0.7` the J = 3/2 domain is 0.26 ≤ T ≤ 0.58, so the
  CLI's own default `--pzz 0.6` is outside it by 0.02 and the run is refused,
  while the same 0.6 is inside the wider spin-1 domain 0.1 ≤ P_zz ≤ 1 and runs
  on ⁶Li. **What the mode costs is measured**
  (`docs/open_items/run_2026-09-06/phase_A_numbers.md` §A3, standard
  configuration, `--pzz 0.5` against the ladder): on **⁷Li nothing physical
  moves** — the rank-2 sector is identically zero, both modes honour `--pz`,
  and σ agrees to 1 ulp (1.97 × 10⁻¹⁶ relative) — while on **⁶Li** the rate
  moves −0.0023527 % and A_∥ by −4.27 × 10⁻⁶ of its own statistical error at
  100 000 events. What moves on both is the **recorded alignment**, 0.4 → 0.5
  at J = 3/2 and 0.409403 → 0.5 at J = 1, which is the divisor of every tensor
  estimator (−20 % resp. −18.1195 % on δ(A_zz) and δ(cos 2φ) at fixed N).
* **At J = 3/2 the (P_z, T) domain is smaller than the spin-1 one.** The four
  populations are fixed uniquely by (1, P_z, T, R₃), so it is a domain and not
  a solver failure: |0.9 P_z + 0.1 R₃| ≤ (1 + T)/2 and
  |0.3 (P_z − R₃)| ≤ (1 − T)/2, i.e. at R₃ = 0
  **1.8|P_z| − 1 ≤ T ≤ 1 − 0.6|P_z|** (so |P_z| ≤ 5/6). The refusal prints
  those edges.

**Why no b₁ is shipped, in one line.** The α–t convolution has been worked out
and measured — one L = 1 partial wave, purely orbital, **2.99 ± 0.02 ×** ⁶Li's
orbital term — but its **sign flips** when the unpolarised backend is moved
from `ToyF2` to MSTW2008 LO, which is the backend the A = 2 gate tells you to
quote. Shipping it would publish a tensor asymmetry whose direction is a flag.

**Two qualifications on the sentence above, because it quotes a research note
and not this tree** (`docs/OPEN_ITEMS_SOLUTIONS.md` §15.4 states them at
length, and this line carried neither until 2026-09-05):

* **The 2.99 is not in the code and was not re-measured here.** It comes from
  `docs/open_items/run_2026-09-03/phase_D_li7_rank2.md` §§5–6; nothing in
  `src/` or `python/` computes a ⁷Li b₁, and `meta["rank2_input"]` says so on
  every ⁷Li run.
* **The Q(⁷Li) gate IS COMMITTED as of 2026-09-06 — and it gates the WAVE
  FUNCTION, not b₁.** `li7_alpha_t_quadrupole` (`b1_nuclear.hpp`) computes
  Q(⁷Li) = **−3.485059 fm²** from `li7_at3.momentum` and reports the ratio
  against **`LI7_QUADRUPOLE_FM2` = −4.06 fm²**, now in `rc.hpp` beside
  `LI6_QUADRUPOLE_FM2` and sourced from the same TUNL A = 5, 6, 7 evaluation
  (Tilley *et al.*, NPA 708 (2002) 3, `Q = −40.6 ± 0.8 mb`) the ⁶Li constant
  comes from: **ratio 0.858389, 14 % low** (0.871265, 13 % low, against the
  −4.00(3) fm² the research note had quoted from memory, which is one of the
  eight determinations N. J. Stone lists without recommending one). Pinned by
  doctest T13 and pytest G8; measured in
  `docs/open_items/run_2026-09-06/phase_B_numbers.md` §B3; author decision
  **D11** is thereby paid. **It licenses no b₁**: what it validates is the α–t
  wave function's ⟨r²⟩ and P-wave character, ⁷Li's rank-2 sector is still
  exactly zero, and the pytest that runs the gate re-asserts that zero in the
  same test. There is also no A = 3 analogue of the A = 2 gate at all: ³H/³He
  are J = ½ and carry no rank-2 structure function, so nothing about the ⁷Li
  convolution can be validated on a lighter system.

The full construction, its numbers and the author decisions it needs are
`docs/OPEN_ITEMS_SOLUTIONS.md` §15 and
`docs/open_items/run_2026-09-03/phase_D_li7_rank2.md`.

## 3. Tagged (⁶Li α, ⁷Li α, d control)

```cpp
PipelineConfig cfg;
cfg.channel  = PipelineChannel::TaggedLi6Alpha;   // or TaggedLi7Alpha / TaggedDeuteronP
cfg.isotope  = channel_isotope(cfg.channel);      // "6Li" / "7Li" / "d"
cfg.beam_config = 1;                              // 10 GeV e x 99.5 GeV/u
cfg.n_events = 400000;
cfg.optics_choice = OpticsChoice::Tagging;        // route label; see below
Pipeline p(cfg, tensor_thirds_plan(0.7, 0.6));

p.for_each([&](const Event& ev) {
  const Particle* alpha = ev.find(Role::Spectator);      // on shell, status 1
  const Particle* dstar = ev.find(Role::StruckCluster);  // OFF shell, status 3
  const Particle* pn    = ev.find(Role::StruckNucleon);  // T1: OFF shell, status 3
  if (rp_tagged(ev, p.optics(), p.pot_config()))         // Roman-Pot mask
    printf("k=%.3f cos=%+.3f m_S=%+g alpha_s=%.3f pT_s=%.3f  N=%d pol=%+g\n",
           ev.kin.k, ev.kin.cos_theta_k, ev.spin.m_struck,
           ev.kin.alpha_s, ev.kin.pt_s, pn->pdg, pn->pol);
});
```

The DIS side is the struck cluster's, not the ion's: `InclusiveKinematicsSource`
wraps an `InclusiveSampler` built on `channel.dis_target` (embedded **deuteron**
for the ⁶Li α tag, quasi-free **triton** for the ⁷Li α tag, free **neutron** for
the d control) and conditions the (x, Q², φ) draw on the struck-cluster
projection m_S through a pure spin category of spin S_c — the C++ shape of
`polligen.tagged.TaggedSampler._pure_category`.

**No inclusive b₁ by default** (`StruckClusterOptions::inclusive_b1 = false`).
In the impulse approximation the embedded deuteron's b₁ *is* the k-integral of
the m-dependent spectator density the sampler already draws from, so an
inclusive b₁ double-counts the tagged tensor asymmetry
(`money_tagged_azz.py`, 2026-08-29). Set `cfg.struck.inclusive_b1 = true` only
for the k-integrated **rate** identity
`Azz(σ₊₁, σ₀, σ₋₁) = tensor_dilution × ⟨Azz⟩_σ`.

### Cluster radial forms: Hulthen (default) or the ANL VMC tables

```cpp
cfg.cluster_wave = ClusterWaveSource::VmcAV18;   // default: Hulthen
```
```bash
python -m lipolgen.cli --channel tagged-alpha --cluster-wave vmc --events 400000
```
```python
cfg = lipolgen.make_config(channel="tagged-alpha", events=400000,
                           cluster_wave="vmc")     # or "hulthen" (default)
```

`Hulthen` is the two-parameter analytic family (`--cluster-beta`, `--p-d`) and
stays the default, so every published number is unchanged bit-for-bit.
`VmcAV18` swaps in the tabulated ANL variational-Monte-Carlo cluster wave
functions for the two **lithium alpha tags** — magnitudes from
`data/vmc/momenta/`, the S–D relative sign from `data/vmc/li6_alpha_d/li6.ad`
— and then **ignores `cluster_beta` and `p_d`**: the shape is the table's, and
P_D(⁶Li) is a property of the wave function (1.935 %, against the 8.67 %
scenario placeholder).  Data
files are found at `$LIPOLGEN_DATA_DIR`, else the compiled-in
`${CMAKE_SOURCE_DIR}/data`; the tables are zero past 5 fm⁻¹ = 0.9866 GeV, so
no spectator is drawn beyond that.

**Since 2026-09-04 the flag also reaches the deuteron control channel** (open
item C5.4).  It used to be **silently ignored** there — the comment said "there
is no VMC d → p+n cluster table, the deuteron IS the cluster", which is true of
the `momenta/` files and false of `data/vmc/deuteron/fdeut.av18`, whose u(k)
and w(k) *are* the p–n relative S and D waves.  `--cluster-wave vmc` now builds
the control from that table at its own P_D = **0.057600**
(`deuteron_av18_p_d()`, +28 % over the scenario `P_D_DEUTERON` = 0.045):

| deuteron control | P_D | `vector_dilution` | `tensor_dilution` |
|---|---|---|---|
| Hulthén (default, bit for bit) | 0.045 | 0.932496109 | 0.959488878 |
| AV18 `fdeut` | 0.0575998920 | 0.913595979 | 0.948146339 |
| change | +28.00 % | **−2.03 %** | **−1.18 %** |

*(The four dilutions were re-measured on the fixed library 2026-09-06 and moved
in the **6th–7th decimal, ≤ 1.4e−6 absolute** — 0.932494769 / 0.959488074 /
0.913594777 / 0.948145618 before, i.e. the two **vector** dilutions move by
1.34e−6 and 1.20e−6 (6th decimal) and the two **tensor** ones by 8.0e−7 and
7.2e−7 (7th).  They are angle-integrated, so `∫Θ₀Θ₂ dc = 0` by L-orthogonality kills
the S–D cross term the sign fix flips and only the 96-cell midpoint-quadrature
residual of that zero survives: measured 0.932496109312 / 0.959488878164 /
0.913595978560 / 0.948146339359, the two percentages **−2.026832 %** and
**−1.182144 %** (−2.026820 % / −1.182136 % before: at five decimals the
**vector** leg moves in the last place, −2.02682 → −2.02683, and the tensor
leg not at all; at the three digits row 9 prints, −2.03 % / −1.18 %, neither
moves).  P_D does not move at all — it is
a norm.  `tests/test_tagged.cpp` T25 pins all four at rtol 1e−6.)*

The relative S–D sign does **not** flip: CDKS fix φ₂ = −W with φ_L = i^L ψ_L,
so the physical deuteron has ψ₂ = +W > 0 at low k, which is what the
positive-definite Hulthén forms already assume for the **stored** ψ — the
*opposite* of the ⁶Li α–d case.  **What that does not license (corrected
2026-09-06):** ψ₂ = +W is not what the partial-wave amplitude sums.
`TaggedModel::build_amp2` consumes **φ_L = i^L ψ_L**, so it needs
φ₂ = i²ψ₂ = **−W**, and until 2026-09-06 it applied no phase at all — the
tagged sector's S–D interference sign was therefore inverted, against
Cosyn–Weiss II Eq. (6.12) (`src/core/tagged.cpp` `build_amp2`, one
`(-1)^floor(L/2)`; `docs/benchmarking/07_cw_sign_investigation.md`).  `rad_[2]`
and `Wave::psi()` are **not** negated, so ψ₂ = +W remains the stored convention
this paragraph is about.  Hulthén stays the default because the control's job is to be the
Cosyn–Weiss tagged limit of the same analytic family the ⁶Li channel is built
from, and because `fdeut.av18` prints no MC errors (so the AV18 row carries no
band of its own).

That rationale — *switching only one thing would make the control and the
channel it controls two different wave-function families* — is what §C5.5b
then had to apply to the opt-in path as well: until 2026-09-04
`--cluster-wave vmc` made the **control's** deuteron AV18 and left the ⁶Li α
channel's **embedded** deuteron on the 0.045 scenario, which is the same split
one level down.  It does not any more; see the paragraph below.

**The ANL Monte Carlo band, `vmc_mc_sigma` / `cluster_vmc_mc_sigma`** (open
item C5.2).  The ANL momentum files print a 1σ error column per point; until
2026-09-04 the readers parsed it and threw it away, so `VmcAV18` had no error
band at all.  `VmcRadial` now carries `dpsi()`, `shifted_by_sigma(n)` and
`norm2_error()` (which returns the **correlated** and the **quadrature** limit,
because one variational walk fixes neither), and the run-level knob is

```cpp
cfg.cluster_wave = ClusterWaveSource::VmcAV18;
cfg.cluster_vmc_mc_sigma = +1.0;    // or -1.0; 0 (default) is bit for bit
```
```python
ch = lipolgen._lipolgen.li6_alpha_channel(
    source=lipolgen._lipolgen.ClusterWaveSource.VmcAV18, vmc_mc_sigma=1.0)
```

It shifts every point of both waves by the same n·σ (fully correlated — the
conservative envelope) and moves P_D by the ratio the shifted norms imply.
`validate()` **refuses** it off a lithium α-tag channel or on `Hulthen`, by the
same rule the b₁ band scales are refused under: a knob that did not run may not
be recorded as if it had.  Measured, it is small: **0.47 % / 1.62 %**
(correlated) on the S/D norms, **1.1 %** on P_D per σ, and **0.02 %** on the
tagged `tensor_dilution` — against 6.6 % for the Hulthén → VMC choice itself.
There is no CLI flag; set the config field.

**The embedded deuteron follows the flag too — since 2026-09-04** (open item
C5.5b).  A ⁶Li α-tag run reads the embedded deuteron's own wave function in
**two** places besides the α–d relative motion: `TaggedChannel::dis_target`,
which is the struck cluster's g₁, and the T1 `BreakupOptions`, which is the
struck-nucleon **momentum and** spin draw — `ClusterBreakup` samples k from
that same deuteron (`src/core/breakup.cpp` `sample_kc_one`).  Until 2026-09-04
**neither** followed
`--cluster-wave vmc`: the relative motion came from the ANL VMC AV18+UX
overlap and the embedded deuteron stayed on the scenario Hulthén
P_D = 0.045, so every polarized tagged-α observable came out **+2.069 %**
high against the AV18 deuteron (P_D = 0.057600, `deuteron_av18_p_d()`) that
belongs to that overlap — two deuteron wave-function families inside one run.
Both now follow it (`DEUTERON_AV18()`; `BreakupOptions::source`).

| ⁶Li α tag, `--cluster-wave vmc` | embedded-deuteron dilution | whole-nucleus reading |
|---|---|---|
| before 2026-09-04 | 0.932500 (scenario Hulthén) | 0.905427 (α–d VMC × scenario d) |
| **now** | **0.913600** (AV18 `fdeut`) | **0.887076** (both AV18) |
| change | **−2.027 %** | −2.027 % |

The Hulthén default is untouched, bit for bit.  ⁷Li deliberately does **not**
move: there is no AV18 A = 3 wave function in this tree, so on the α–t channel
the flag selects the relative motion alone.

**The cost is not only on the polarized observables.**  The same deuteron
supplies the T1 struck-nucleon MOMENTUM draw, so the `struck_virtuality`
column — unpolarised and plan-independent — moves on this flag as well.
Measured on `--cluster-wave vmc --events 400 --seed 4242` (tagged-6Li-alpha,
`tensor-thirds`) against the same command on the pre-fix build: **378 of 400
events** differ, ⟨p² − M²⟩ **−0.0567 → −0.0599 GeV²**, σ ×1.54; `k`, `weight`
and every other T0 column are bit-identical, and the Hulthén default is
untouched there too.

**And the inclusive constants do NOT follow the flag.** Under
`--cluster-wave vmc` a tagged row and an inclusive row of the same programme
describe ⁶Li with wave functions whose vector dilutions differ by **+11.61 %**
and whose rank-2 transfers differ by **+6.58 %** (open item C5.5).  Since
2026-09-04 `validate()` **refuses** `cluster_wave` on `--channel inclusive` and
`--channel coherent`, where it is never read: accepting it unread was how that
11.61 % reached someone who had asked for a "VMC ⁶Li" row and got the shipped
0.811228 with no warning.  Substituting the VMC P_D into
`LI6_CLUSTER_POLARIZATION` is **not** the fix: it gives 0.905427, 6.8 % *above*
the ab-initio six-body VMC 0.848, where the shipped 0.811228 is 4.3 % below it.
**Quote no inclusive ⁶Li polarization without the band 0.81 … 0.91.**

What it changes (40 k events, 10 × 99.5 GeV/u; full tables and the derivation
in `docs/open_items/vmc_reconciliation.md`):

| | Hulthen β = 0.30 | VMC AV18 |
|---|---|---|
| ⁶Li ⟨k⟩ / P(k>0.2) / P(k>0.45) | 0.1219 / 0.1476 / 0.0157 | 0.1225 / 0.2410 / 0.0019 |
| ⁷Li ⟨k⟩ / P(k>0.3) / P(k>0.45) | 0.2893 / 0.3527 / 0.1555 | 0.1864 / 0.2111 / 0.0114 |
| ⁶Li P_D | 0.0867 | 0.01935 |
| ⁶Li tag fraction, YR high-acceptance | 0.0264 | 0.0348 |
| ⁶Li tag fraction, tagging optics | 0.2551 | 0.2486 |
| ⁷Li tag fraction, YR high-acceptance | 0.9730 | 0.9981 |
| ⁶Li A_zz^tag at k = 0.20 GeV | **−1.207** | **−0.519** |

> **A_zz^tag row corrected 2026-09-06.**  It read `+0.845 | +0.452` until the
> S–D interference sign of the tagged amplitude was fixed
> (`src/core/tagged.cpp` `build_amp2`, one `(-1)^floor(L/2)`: the partial-wave
> sum needs φ_L = i^L ψ_L and had no phase at all).  Measured at k = 0.1979 GeV
> — the grid cell nearest 0.20 — on the acceptance-weighted curve at
> `default_configs("6Li")[1]` and `yr_optics(..., high_acceptance = true)`,
> `n_phi = 32`: **−1.2069** (Hulthén β = 0.30) and **−0.5191** (VMC AV18),
> against +0.8450 and +0.4518 before.  The magnitude of the change is the
> whole observable: A_zz^wf lives in [−2, 1] and the sign flip moves it by up
> to 2.74.  ⟨k⟩, P(k>·) and P_D are spin-blind and unchanged to the digits
> shown, and the spin-blind accepted rate is **unmoved across the fix to
> ≤ 4.0e−16 relative (1–2 ulp; the VMC channel is bit-identical)** — measured
> as the acceptance-weighted Σ_M n_M k² at `n_phi = 32`, giving an accepted
> fraction of **0.0246759321488** (Hulthén β = 0.30) and **0.0338102276258**
> (VMC AV18). Digits beyond those are summation-order dependent, not physics.
>
> **The two ⁶Li tag-fraction rows DID move, and the claim that they did not is
> withdrawn (2026-09-06, verification pass).**  They read `0.0249` and `0.2530`
> until this edit.  They are 40 k-event SAMPLES at the `tensor-thirds` fill —
> a category average whose EXPECTATION is the uniform-M mix integral and IS
> invariant, measured unmoved to all 15 digits (**0.024675932148828** Hulthén
> β = 0.30, **0.033810227625842** VMC AV18, the equal-thirds average
> reproducing it to +3.1e−6 (a per-M norm-residual effect — Σ n_M k² differs by 5e−5 between M = 0 and ±1 — not summation order; the pooled category average equals the uniform mix to 2.2e−16 on both builds)) — but whose SAMPLE is re-drawn, ~77 % of events
> taking a different (k, cos θ_k); at 40 k events the two readings sit 1.4σ
> (YR high-acceptance, 0.0249 → 0.0264) and 0.7σ (tagging, 0.2530 → 0.2551)
> apart, on binomial σ_diff of 1.1e−3 and 3.1e−3 for two independent samples.  A single tensor-polarised fill is not invariant at all: measured on
> the ⁶Li Hulthén channel at these optics, the `tensor-thirds` categories go
> **0.028127 → 0.018403** (azz±, −34.6 %) and **0.017775 → 0.037222** (azz0,
> +109 %), and the `--pz 0.7` max-entropy ladder of `helicity-flip`
> (0.751567 / 0.196866 / 0.051567) goes **0.027030 → 0.020396**, −24.5 %.
> ⁷Li does not move at all.
> Reason and derivations:
> `docs/benchmarking/07_cw_sign_investigation.md`; before/after tables:
> `docs/open_items/run_2026-09-06/phase_CW_numbers.md`.

**Roman-Pot tag fractions, Hulthén β = 0.30 vs VMC, at the three configurations
(40000 events, seed 20260829):** `validation/vmc_tag_fractions.py --events
40000 --configs 0,1,2` (5×41-class, 10×100-class, 18×275-class nominal design
points; 18×275 rigidity-caps ⁷Li's ion energy below ⁶Li's, so it is quoted
separately per isotope).

| isotope | optics | 5 × 40.8 GeV/u Hulthén | VMC | 10 × 99.5 GeV/u Hulthén | VMC | 18 × 137.5/117.9 GeV/u Hulthén | VMC |
|---|---|---|---|---|---|---|---|
| ⁶Li | YR high-acceptance | 0.0301 | 0.0365 | 0.0264 | 0.0348 | 0.0279 | 0.0348 |
| ⁶Li | tagging optics | 0.3451 | 0.3185 | 0.2551 | 0.2486 | 0.3145 | 0.2911 |
| ⁷Li | YR high-acceptance | 0.9660 | 0.9981 | 0.9730 | 0.9981 | 0.9787 | 0.9981 |
| ⁷Li | tagging optics | 0.9805 | 0.9993 | 0.9927 | 0.9992 | 0.9941 | 0.9992 |

**This table was REGENERATED on the fixed library, 2026-09-06** (the command
above, unchanged, on the post-`build_amp2`-fix build).  The ⁶Li rows read
`0.0286 | 0.0365 | 0.0249 | 0.0348 | 0.0266 | 0.0349` and
`0.3410 | 0.3175 | 0.2530 | 0.2485 | 0.3115 | 0.2905` before, and a scratch
build of the pre-fix `HEAD` reproduced the six Hulthén cells exactly when the rows were checked (two VMC tagging cells differ in the last printed digit: 0.2486 vs 0.2485 at 10 × 99.5, 0.2906 vs 0.2905 at 18 × 275)
on 2026-09-06; **every ⁷Li row is unchanged, digit for
digit**, because ⁷Li is one L = 1 wave and does not move at all.  The ⁶Li rows
are a `tensor-thirds` category average, so what moved is the SAMPLE and not the
expectation — the uniform-M mix accepted fraction is unmoved to all 15 digits
(0.024675932148828 Hulthén β = 0.30, 0.033810227625842 VMC AV18) — but the
published digits are no longer the pre-fix ones, which is why they are
regenerated rather than annotated.  For what a single tensor-polarised fill
does (up to −24.5 % at the CLI's own default fill), see the correction box
above and `docs/open_items/vmc_reconciliation.md`.

The cluster-wave systematic is quoted as the difference between a
`--cluster-wave hulthen` run and a `--cluster-wave vmc` run; the Hulthén β
band is retired as a systematic (the `--cluster-beta` knob stays).

Reach the tables directly when you want them: `VmcRadial`,
`vmc_from_overlap_k`, `vmc_from_overlap_r` (a Fourier–Bessel transform of an
r-space overlap, for the r-only tables) and `vmc_from_momentum` are all bound
in Python, and `TaggedModel::radial_table(l)` gives the normalized
ψ̂_L(k)·√P_L the sampler actually draws from.

### Tier T1 — the struck cluster is resolved (the default)

`PipelineConfig::tier` is `Tier::T1`, so the tagged record does not stop at
the cluster. `ClusterBreakup` (`breakup.hpp`) draws the internal relative
momentum from the cluster's own wave function and writes

* one `Role::StruckNucleon` (status 3, **off shell**, `pdg` 2212/2112, `pol`
  a sampled ±1 helicity label), and
* its `Role::PartnerSpectator` fragments (status 1, **on shell** at their
  AME2020 masses, 10-digit ion codes for nuclei).

| channel | struck cluster | partners | model |
|---|---|---|---|
| ⁶Li α | embedded d | 1 nucleon | S + D Hulthén at `P_D_DEUTERON`, the m_S-dependent \|A_{m_sc}\|²; species `F2p : F2n`; nucleon spins from the pair CG factor, so ⟨P_N⟩ = (1 − 3/2 P_D) m_S |
| ⁷Li α | quasi-free t | d, or two n | sequential two-body at the AME2020 S_n(³H) = 6.2572 MeV / S(p+n+n) = 8.4818 MeV, nn split at 1/\|a_nn\| = 10.44 MeV — **crude and flagged** (plans/05 5.D) |
| d control | already a nucleon | none | relabelled only |

The impulse-approximation rule is applied twice: **the spectators are
physical, the struck object is not.** Partners go on shell, the struck
nucleon takes `P_X − Σ p_partner`, and no fragment is ever recoil-corrected.
Mean virtuality `p_N² − M_N²`: −0.057 / −0.124 / −0.035 GeV².

Balance at T1: **whole nucleus**,
`k + P_ion = k' + p_spec + Σ p_partner + X`, with `X = k + p_N,struck − k'`
the per-nucleon remainder. The tagged spectator is on shell at its AME2020
mass and the struck cluster `P_X = P_ion − p_spec` stays on the record as
documentation. Charge closes.

Set `cfg.tier = Tier::T0` for the pre-T1 record
(`k + P_ion = k' + p_spec + X`, `X = k + P_X − k'`, no partners, no struck
nucleon), which is what every tagged number published before this tier was
made with — every T0 quantity is bit-identical between the two apart from the
0.02 % of draws T1's own timelike-X rejection redraws.

Naming the struck nucleon is also what makes the **T2 chain conserve with no
caller-side hook** (§6, `docs/T2_CHAIN.md` §1a).

### Triton spectral function: sequential Hulthén (default) or Ciofi–Simula

```cpp
cfg.triton_sf = TritonSfChoice::CiofiSimula;    // default: Hulthen
```
```bash
python -m lipolgen.cli --channel tagged-7Li-alpha --triton-sf ciofi-simula --events 400000
```
```python
cfg = lipolgen.make_config(channel="tagged-7Li-alpha", events=400000,
                           triton_sf="ciofi-simula")  # or "hulthen" (default)
```

`hulthen` keeps the sequential two-body triton decay of the table above
bit-for-bit.  `ciofi-simula` replaces it, on the **⁷Li α tag only**, with the
Ciofi degli Atti–Simula spectral function (`triton_sf.hpp`): the struck
neutron's momentum is drawn from n₀(k) + n₁(k) and the **branching between a
bound deuteron remnant and a (p n) continuum remnant is the k-dependent ratio
n₀/(n₀ + n₁)** — never an assumed constant — whose k-integral is the
³He(e,e′p)d spectroscopic factor S₀ = 0.6525.  Three channels instead of two:

| channel | weight | remnant | measured fraction |
|---|---|---|---|
| struck n → n + d | n₀/(n₀+n₁) | bound d, E = 0 | 65.3 % of struck n |
| struck n → n + (pn) | n₁/(n₀+n₁) | continuum at the pn ¹S₀ pole, 8.31 MeV | 34.7 % of struck n — **the channel the sequential model does not have** |
| struck p → p + (nn) | always | continuum at the nn pole, 10.44 MeV | every struck p |

⟨k⟩ moves from 133 MeV (n + d) / 145 MeV (p + nn) sequential to 102 MeV on
the bound channel and 126 MeV over all struck neutrons.  The species stays
`Z F2p : N F2n`, the impulse-approximation rule is unchanged (partners on
shell, struck nucleon absorbs the difference — conservation is untouched:
whole-record closure < 1e-9 through T1 **and** the T2 PYTHIA chain, both
options, `tests/test_triton_sf.cpp` / `tests/test_t2.cpp`), and `sample()`
consumes a fixed six uniforms on every branch so the stream stays aligned
whichever channel comes out.  The `Pipeline` builds the `CiofiSimulaTriton`
itself at the run's own `--cluster-beta`; a C++ caller who sets
`cfg.breakup.triton_sf` directly (e.g. a `CiofiSimulaOptions` with
`n1_scale`, `proton_n1_only`, or a future Faddeev table behind
`TritonSpectralFunction`) keeps their own object.  Coefficient provenance —
CS PRC 53 (1996) 1689, Eq. (74)/(76), Tables A.1/A.3 — and the BeAGLE
n₀-only caveat: `docs/CONVENTIONS.md` and the `triton_sf.hpp` header.

### FSI of the DIS debris with the spectator — a weight, never a shift

```cpp
cfg.fsi = PipelineFsi::GlauberCluster;     // default: Off = today's PWIA
cfg.fsi_sigma_mb = 40.0;                   // band 20–40 mb; run BOTH ends
```
```bash
python -m lipolgen.cli --channel tagged-alpha --fsi glauber-cluster \
                       --fsi-sigma-mb 40 --events 400000
```
```python
cfg = lipolgen.make_config(channel="tagged-alpha", events=400000,
                           fsi="glauber-cluster", fsi_sigma_mb=40.0)
```

The eikonal (Cosyn–Weiss/Glauber) rescattering of the hadronic debris X on
the tagged cluster, applied as a **multiplicative `Event::weight`** — the
FSI/IA density ratio at the drawn (k, cos θ_k), built once at setup
(`GlauberFsiWeight`, `fsi.hpp`) on a (k_z, k_T) grid — the kernel transfers
transverse momentum only, so the table is smooth and even in k_z there; a
spec-literal (k, cos θ_k) grid would tabulate the same function on skewed
axes — and read per event by bilinear interpolation.  Costs ~0.2 s at setup
and nothing measurable per event; no RNG is consumed, so an FSI-on run is
**bit-identical to the FSI-off run in every four-vector** and only the
weight column moves (`tests/test_fsi.cpp` proves this event by event).

Three rules, from the header, in order of importance:

* **A weight, never a shift.** The spectator four-vector is the measurement;
  moving it would break the "never recoil-correct the light spectator" rule
  and the whole-nucleus balance. `Off` is the plane-wave impulse
  approximation bit for bit.
* **Spin independent by construction.** The weight is the m-summed
  (unpolarized) shape distortion; nothing constrains the spin dependence of
  the rescattering (Cosyn–Weiss VI C, open question). Quote it as an
  **unpolarized-shape systematic, never as a correction to A_zz**.
* **Band σ_XN over 20–40 mb; never quote one row alone.** 40 mb is the free
  hadron; at EIC formation lengths the 20 mb row is arguably the realistic
  one. `--fsi-sigma-mb` sets the row; run both ends as the systematic.

The run summary logs the survival probability (∫w dΓ/∫dΓ, what
`weight_normalised` divides by — the distortion is ~73 % absorptive, so it is
well below 1 and that is physical): 0.520 at 40 mb, 0.671 at 20 mb on the
⁶Li α tag (S+D). `glauber-cluster` (the default variant) shadows the X–α
cross section over the α's own Gaussian profile — σ_Xα = 131.0 mb at
σ_XN = 40, B_α = 27.2 GeV⁻², σ_el/σ_tot = 0.269 vs the measured α-p 0.258 —
while `glauber-nucleon` is the **unshadowed A·σ_XN = 160 mb single-scattering
limit**: more absorptive point by point at low k, though its integrated
survival lands *above* the cluster variant's (the unshadowed quadratic gain
term feeds strength back into the tag; `fsi.hpp` header, pinned in
`tests/test_fsi.cpp`). The weight reaches HepMC3 `weights()[0]` and the npz
`weight` column with no further wiring; `Pipeline::fsi_weight()` exposes the
model (`survival()`, `sigma_eff_mb`, the profile) for printing.

**The two variants are NOT one** — measured 2026-09-04 (D3), because the
open-items inventory said they were. On ONE event stream reweighted three
ways (`--channel tagged-6Li-alpha --events 20000 --seed 1234`, plan
`tensor-thirds` at P_z = 0.7 / P_zz = 0.6 / P_e = 0.7, at the default
σ_XN = 40 mb — one stream at one end of the mandatory 20–40 mb band; the
columns `k`, `cos_theta_k`, `phi_k`, `x`, `q2` bit-identical across the three
files):

| | `off` | `glauber-cluster` | `glauber-nucleon` |
|---|---|---|---|
| Σ weight (20 000 events) | 20000.0 | 10419.07 | 11632.10 |
| Σw/Σw_off (the sample mean weight) | 1 | 0.520954 | 0.581605 |
| `survival()` — the model's own grid-integrated one, `meta["fsi_survival"]` | 1 | 0.520239 | 0.582899 |
| per-event weight, min … max | 1 | 0.0208 … 1.4043 | 0.2139 … 6.9523 |

Per event, w_nucleon/w_cluster has median 0.897 and reaches **68.52**, and
**99.50 % of the events differ by more than 1 %**. Different by construction,
not by parameter choice: the algebraic identity that *does* hold (the Ciofi
degli Atti–Kaptari per-nucleon product on an uncorrelated density collapses
onto the cluster form) is about a code path this library does not have, which
is exactly why the shipped `glauber-nucleon` is the single-scattering limit
instead.

**What the run records.** Since 2026-09-04 an FSI run's npz `meta` carries
`fsi`, `fsi_sigma_mb`, `fsi_sigma_cluster_mb`, `fsi_sigma_cluster_el_mb`,
`fsi_survival`, `fsi_clipped_grid_fraction`, `fsi_formation_ramp` (and the
four ramp anchors when the ramp ran). Before that the three runs above
produced **identical `meta` dicts** while their total rate differed by 48 % /
42 % — which made the mandatory 20–40 mb band unquotable, because the two ends
of it were indistinguishable files. An `--fsi off` npz carries exactly the old
key set: the block is conditional.

## 4. Coherent ⁶Li

```cpp
PipelineConfig cfg;
cfg.channel = PipelineChannel::CoherentLi6;
cfg.isotope = "6Li";
cfg.n_events = 300000;
cfg.coherent.f0      = 0.04;     // band 0.02 - 0.08     (SCENARIO)
cfg.coherent.slope_b = 50.0;     // band 40 - 60 GeV^-2  (SCENARIO)
cfg.coherent.amp     = 0.01;     // flat cos2phi at P_zz = 1, band 3e-3 - 1e-2
cfg.coherent_t_max   = COHERENT_T_MAX_DEFAULT;   // 0.2; larger THROWS (below)
cfg.coherent_xpom.m_x_min   = 1.2;   // smallest diffractive mass [GeV]
cfg.coherent_xpom.x_pom_max = 0.1;   // upper edge of the diffractive region
cfg.optics_choice = OpticsChoice::Tagging;   // the YR envelope tags nothing
Pipeline p(cfg, tensor_thirds_plan(0.7, 0.6));

p.for_each([&](const Event& ev) {
  const Particle* recoil = ev.find(Role::IntactRecoil);
  const double phi_t = std::atan2(recoil->p.py, recoil->p.px);
  // phi_t is drawn flat and the modulation rides on ev.weight, so
  // 2 <w cos 2(phi_t - phi_S)> estimates the cos 2phi coefficient
});
```

The coherent yield is `f_coh(x) × the UNPOLARIZED inclusive rate`
(`coherent.project_coherent`: `n_coh = n_events × f_coh(x)`), so σ is spin
independent and the tensor signal lives entirely in the **recoil azimuth**,
`1 + c₂ cos 2(φ_t − φ_S)`. Folding the inclusive `w_avg` in on top would count
the polarization twice. Balance: `k + P_ion = k' + P_recoil + X`, with the
neutral diffractive system X carrying the remainder.

**The pomeron.** `CoherentXpomModel` draws the per-nucleon fraction
`x_P = (M_X² + Q²)/(W² + Q²)` log-uniformly on `[x_P(M_X,min), x_pom_max]`, so
every event carries a diffractive mass of at least `m_x_min`
(`COHERENT_MX_MIN_DEFAULT` = 1.2 GeV, above the ρ/ω/φ region — the exclusive
vector-meson channel is a different process, deliberately not generated;
`coherent.hpp` documents the choice and the hadronization veto table behind
it). The recoil is then **solved**, not approximated:
`(k + P_ion − k′ − P_recoil)² = M_X²` with `P_recoil² = M_A²` and
`p_T = √|t|` is a quadratic in the recoil's light-cone plus momentum, so the
balance closes to rounding and X is timelike by construction. The record
carries `kin.x_pom`, `kin.beta_pom = x/x_P` and `kin.m_x2`; the nucleus loses
`x_P/A` of its light-cone momentum, which keeps the recoil rigidity in
[0.979, 1.000] — inside the near-beam band, where it has to be to be tagged at
all. Cells that cannot fit `m_x_min` below `x_pom_max` carry no coherent rate.

**One thing throws at setup rather than biting later.** `coherent_t_max` beyond
the range where `1 + c₂ cos 2(φ_t − φ_S)` stays positive (0.245 GeV² at
P_zz = −2, 0.495 at P_zz = +1) is refused — `CoherentScenario::positivity_margin`
is the coherent twin of `InclusiveKernel::positivity_margin`.

### The |t| ceiling — why 0.2, and which reason survives a change of `eps_b0`

```bash
python -m lipolgen.cli --channel coherent --coherent-t-max 0.2 --events 200000
```

`COHERENT_T_MAX_DEFAULT` = 0.2 GeV² **stays**, and since 2026-09-04 (D5) the
reason it is written down is the **anchor range**, not positivity. The two are
not interchangeable, and the header used to present them as agreeing:

* **Primary, and independent of every knob.** The deformation mechanism is
  scaled from Mäntysaari *et al.*'s polarized-deuteron a₂, and
  `mantysaari_a2_deuteron()` carries **four digitised rows, |t| = 0.05, 0.10,
  0.20, 0.30**. The fit is linear in |t| and exact only as |t| → 0, so 0.2 is
  inside the input and 0.5 is outside it. This is a property of the input
  table; it does not move when `eps_b0` moves.
* **Secondary, and contingent on `eps_b0`.** The linear c₂ crosses −1 at
  |t| = 0.245 (P_zz = −2). That number is now **derived**, by
  `CoherentScenario::t_positivity_edge(pzz)`:

```python
sc = lipolgen.CoherentScenario()
sc.t_positivity_edge(-2.0)          # 0.245  -- the shipped eps_b0 = -0.08
sc.eps_b0 = -0.0070                 # the MEASURED 6Li quadrupole row, ROUNDED
sc.t_positivity_edge(-2.0)          # 2.8000 on that rounded input; the DERIVED
                                    # eps_b0 = -0.0070024 gives 2.7990.  Quote
                                    # 2.80 GeV^2, never "2.8000" -- 9.3x
                                    # outside |t| <= 0.30 either way.
```

  So **positivity stops binding the moment `eps_b0` is corrected, and the
  anchor range does not.** A positivity-derived ceiling would license
  extrapolating a linear-in-|t| fit ten times past its data — less honest than
  the fixed number, not more. (`eps_b0` is 11.4× the measured ⁶Li quadrupole;
  see "ΔB, `eps_b0`…" below and STATUS.md decision row 8.)

`positivity_margin` / `check_positivity` are unchanged and still throw: they
are the **guard** on the truncated weight the sampler actually uses, catching
an author who raises `t_max` or `eps_b0` past where it stops being a density.
That is a different job from justifying the constant.

**The ceiling costs no rate.** Measured on 200 000 generated events at the
shipped defaults (seed 99, plan `tensor-thirds` at P_z = 0.7 / P_zz = 0.6, 4
threads): ⟨|t|⟩ = 0.019963 GeV², max |t| = 0.197575, and the
fraction above 0.05 / 0.10 / 0.15 / 0.20 is 0.082015 / 0.006550 / 0.000455 /
0. `sample_t` **renormalises** on `[0, t_max]`, so the truncation does not
lose rate — it redistributes exp(−B·t_max) = **4.5e−5** of it. The ceiling is
a statement about where the model is defined, not a rate cut.

`--coherent-t-max` (new, 2026-09-04) reaches it from the command line, and a
coherent run's npz `meta` now records `coherent_t_max`, `coherent_slope_b`,
`coherent_eps_b0`, `coherent_amp`, `coherent_f0`, `coherent_m_x_min`,
`coherent_x_pom_max`, `coherent_weighted_azimuth` and the derived
`coherent_t_positivity_edge_pzz_m2`. Before that, a caller who moved
`coherent_t_max` changed the entire |t| spectrum, the tag acceptance and every
c₂ in the file, and nothing in the file said so.

### ΔB, `eps_b0`, and what the shipped scenario assumes about Q(⁶Li)

**ΔB is defined once**, at the `eps_b0` declaration (open item O4, closed
2026-09-04; before that the header used the symbol in three docstrings and
defined it nowhere):

    |F_m(|t|, Φ)|² = exp(−|t| [B + ΔB_m cos 2(Φ − Φ_S)]),   ΔB_m = δ_m/2

with δ_m = ⟨x²⟩ − ⟨y²⟩ per nucleon in state m.  Expanding to first order
against the anchor's `1 + 2 a₂ cos 2Φ` normalisation gives
**a₂(m) = −(ΔB_m/2)|t| = −(δ_m/4)|t|**, which is `a2_m_state` exactly, so

    eps_b0 ≡ δ_{±1}/B = +2 ΔB_{±1}/B = −1 × ΔB₀/B.

The pre-2026-09-04 label "relative slope modulation ΔB₀/B of the m = 0 state"
was off **by a sign**, not by a factor 2; the code was right.
`CoherentScenario::delta_b_m(m)` and `slope_at_azimuth(φ, m)` are the only
places ΔB is computed.

**Ask the scenario what it assumes.**  `quadrupole_from_a2_slope`
(`cluster_config.hpp`) is the exact inverse of `a2_from_quadrupole`:

```python
sc = lipolgen._lipolgen.CoherentScenario()
q_matter = lipolgen.quadrupole_from_a2_slope(sc.a2_m_state(1.0, 1), 6)
0.5 * q_matter            # -0.9345 fm^2 -- the CHARGE quadrupole it implies
```

At the shipped `eps_b0` = −0.08, `slope_b` = 50 that is **11.42×** the measured
`LI6_QUADRUPOLE_FM2` = −0.0818 fm² and 1.52× even the α+d model's −0.615, and
a₂(±1, |t| = 0.3) = **+0.300** — larger in magnitude than the **deuteron's**
own digitized −0.28, for a nucleus whose quadrupole is 3.5× smaller.  The
honest ⁶Li band at this B, from `quadrupole_band_fm2()` through the same map:

| Q_charge assumed | `eps_b0` at B = 50 | a₂(±1, 0.3) | \|t\| at \|c₂\| = 1, P_zz = −2 |
|---|---|---|---|
| measured −0.0818 | **−0.0070** (rounded) | +0.0263 | 2.80 GeV² |
| GFMC −0.20(6) | **−0.0171** (rounded) | +0.0642 | 1.146 |
| α+d model −0.6154 | **−0.0527** (rounded) | +0.1976 | 0.372 |
| **shipped −0.08** | **−0.08** (exact input) | **+0.3000** | **0.245** |

**Read the digit count with the table.** The top three `eps_b0` are ROUNDED and
each row is arithmetic on the rounded value, so the measured-Q row is quotable
as **2.80 GeV²** — three figures, which is all a two-figure input supports —
and **never as "2.8000"**. The DERIVED values
(`ClusterConfigSampler::quadrupole_band_fm2()` through `a2_from_quadrupole` at
B = 50, amp = 0.01; measured 2026-09-05) are `eps_b0` = −0.0070024 / −0.0171207
/ −0.0526846 with P_zz = −2 edges 2.799047 / 1.144810 / 0.372025 — the same
**2.80**, and 2.7990 is the literal `tests/test_coherent.cpp` asserts. This
page printed a four-figure "2.800" with no qualification until 2026-09-05.

The last column is `CoherentScenario.t_positivity_edge(-2.0)`, and its spread
is why `COHERENT_T_MAX_DEFAULT` = 0.2 is **not** justified by positivity: the
edge is a consequence of the oversized `eps_b0`, not of the target, and it
stops binding as soon as `eps_b0` is corrected.  The ceiling's stated reason
is the anchor range instead (see "The |t| ceiling" above).  **The default is
deliberately unchanged** (it is pinned in
`validation/reference/coherent.json`), and the cost is: every *generated*
coherent tensor number is **11.4×** the measured-quadrupole expectation.
**Never publish a single `eps_b0` row** — band it, and say which quadrupole the
row assumes.  Finally, **`eps_b0` and `slope_b` are not independent**: every
observable uses the product δ = `eps_b0`·B, so scanning `slope_b` over {40, 60}
at fixed `eps_b0` moves a₂ by ±20 % for no physical reason.  Band δ.
Numbers and the decision: `docs/open_items/run_2026-09-03/phase_C_numbers.md`
§C4.

**A hadronizer on this channel is an ordinary configuration** (since
2026-08-30; the old refusal and its `hadronize_coherent` opt-in are gone).
The record names its own T2 target — `Role::Pomeron`,
`P_IP = P_ion − P_recoil` — and `PythiaBridge` hadronizes the γ*–Pomeron
system on a PYTHIA Pomeron beam (`Beams:idA = 990`) through the same
surrogate as every other channel, so the whole record conserves exactly
(`docs/PYTHIA_BRIDGE.md` §12). The knobs sit on the bridge, not the pipeline:
`PythiaBridgeOptions::coherent_t2 = {Pomeron, Off}` (Off skips the third
PYTHIA instance and leaves coherent records at T0, still conserving), plus
`pom_set` / `pom_rescale` for the Pomeron PDF; on the command line,
`--coherent-t2 pomeron|off`, `--pom-set`, `--pom-rescale`.

### `--pom-set` — the T2 tier's largest model systematic, MEASURED

```bash
for s in 3 4 5 6 7 8 9 10 12 13 14 15; do
  python -m lipolgen.cli --channel coherent --hadronize --pom-set $s \
                         --events 20000 --seed 4242 --hfs-npz pom_$s.npz
done          # then take the ENVELOPE over the twelve files
```

Scanned 2026-09-04 (D4) over **every set PYTHIA 8.317 ships (1–15), 20 000
coherent events each at ⁶Li config 1, seed 4242** — 17 s wall for the whole
scan, so cost is not a consideration. Set 11 has been refused by the bridge
constructor since (below), so what reproduces today is the other **fourteen**,
in 14.3 s.

**It moves nothing at T0.** The columns `t`, `x_pom`, `q2`, `x` and the event
`weight` are **bit-identical across those fourteen sets** — 1–10 and 12–15,
one md5 `ffd35a3a62b591c547e9ca2ac4301b5d` over the lot, re-measured
2026-09-05; set 11 gave the same md5 in the original scan, but this page said
"all fifteen" until 2026-09-05 and the count that reproduces is fourteen.
The Pomeron PDF enters only the flavour draw and PYTHIA's backward evolution,
while |t|, x_P, M_X and the rate are fixed upstream by `CoherentSampler` /
`CoherentXpomModel`. **A `PomSet` band on M_X, |t|, x_P or σ is identically
zero by construction; do not quote one.** The band is on the **hadronic final
state**:

| observable | set 6 (default) | min | max | band about the default |
|---|---|---|---|---|
| **⟨n_charged⟩** | 3.964 ± 0.016 | 3.860 (set 9) | 4.355 (set 5) | **−2.6 % / +9.9 %** |
| ⟨n_hadrons⟩ | 8.521 | 8.300 (10) | 9.394 (5) | −2.6 % / +10.2 % |
| ⟨p_T⟩ per particle | 0.342 GeV | 0.322 (5) | 0.379 (10) | −5.9 % / +10.9 % |
| **kaon fraction** | 0.0646 | 0.0359 (9) | 0.1002 (10) | **−44 % / +55 %, factor 2.8** |

over the **twelve genuine diffractive-PDF fits** (3–10, 12–15), 20 000
coherent events per set at ⁶Li config 1, seed 4242. Sets 1 (a
Q²-independent toy) and 2 (π⁰ densities) are not Pomeron fits and are outside
the band — **measured, that exclusion is worth 2 points of band**: including
them moves the ⟨n_charged⟩ low edge from −2.6 % to **−4.7 %** (set 2, 3.7773)
and the ⟨n_hadrons⟩ low edge to −7.0 %, leaving the high edges and the ⟨p_T⟩
and kaon bands where they are. So the numbers above are twelve-fit numbers and
must not be quoted against "all 15 sets". **Set 11 is refused** —
`PomHISASD` needs PYTHIA's `setXPom`, which this bridge never calls, and 100.00 % of its events took the charge-democratic
e_q² fallback, so recording `pom_set = 11` would record a knob that did not
run.

Four rules for using it:

* **It is an envelope over re-runs, one npz per set.** The set changes the
  final state event by event; no per-event weight maps one set onto another.
  `meta["pom_set"]` (new, 2026-09-04) is what tells the files apart.
* **Set 6 is the only LO H1 set**, so the band mixes LO and NLO DPDFs used in
  an LO Monte Carlo. Defensible for a systematic envelope, indefensible for a
  central value — say which one a plot is.
* **Set 5 (H1 2007 Jets) makes open charm the default never makes**: 29.98 %
  of its events change when `include_charm` is turned off, against 10.4–11.1 %
  on GKG18 (12–15) and *exactly* 0 on 3, 4, 6, 7, 8, 9, 10 (all three H1 2006
  sets — Fit A NLO, Fit B NLO *and* Fit B LO — carry no charm or bottom at any
  (β, Q²)). Flag it: it is the same signature the `2e404b8` bug faked.
* **20 000 events per set is enough.** Statistical error on ⟨n_charged⟩ is
  0.017 there and the seed-to-seed scatter over 4242 / 777 / 31337 is ≤ 0.057,
  against a 0.50 spread across sets.

**The light-only e_q² flavour fallback is not a systematic.** A default
coherent run takes it on ~20 % of events, and the run banner and
`meta["pom_flavour_fallback_frac"]` now say so. Measured: raising
`q2_pdf_min` from 1.0 to 1.75 drives the share to 0.00 % (sets 6, 3) or a few
per cent (12: 19.45 → 2.65, 13: 20.10 → 3.23, 15: 21.65 → 5.12) — **with one
exception that must travel with the sentence: on set 4 (H1 2006 Fit B NLO) it
stays high, 35.52 % → 29.32 %** — and leaves the final state **bit-identical**
in every one of those six cases, set 4 included. Every Pomeron DPDF carries a
single light-quark singlet, so e_q²·xf_q ∝ e_q² exactly over the light
flavours and the "fallback" *is* the true draw; that is why set 4's stubborn
29 % costs nothing either. (Measured 2026-09-04, re-measured 2026-09-05:
4 000 coherent events per point at ⁶Li config 1, seed 4242, diffing the whole
`pid`/`p4` final state.) (At
`q2_pdf_min` = 3.0 the GKG18 sets do move, but that is a **charm** effect —
the clamp lifts charm above threshold — and is bit-identical again with
`include_charm = false`.) The earlier claim that 1.75 removes the fallback
"at the cost of clamping every flavour weight to that Q²" was true about the
counter and **false about the cost**: it changes nothing.

## 5. Far-forward routing

The route is **not stored** on the event; it is recomputed from the record, so
one generated sample can be priced at several envelopes:

```cpp
const double p_u = p.beam_config().ion_momentum_per_nucleon;
Optics ha  = optics_for(OpticsChoice::YellowReportHighAcceptance, "6Li", p_u);
Optics tag = optics_for(OpticsChoice::Tagging, "6Li", p_u);      // de-squeezed
int   r    = route_of(ev, ha, p.pot_config());   // kRouteRomanPots, kRouteB0, ...
bool  seen = rp_tagged(ev, tag, p.pot_config()); // main window + near-beam tail
```

`OpticsChoice::TaggingLegacyLevers` prices the tagging point with the 18×275 pot
levers everywhere — the way every tagging number published before 2026-08-29 was
made. `OpticsChoice::Custom` uses `PipelineConfig::optics` verbatim.

Measured ⁶Li α tag at 10 × 99.5 GeV/u (400 k events): **0.0248** at the Yellow
Report high-acceptance optics, **0.3046** at the tagging optics with the legacy
levers, **0.2547** with the per-configuration levers — against the Python's
0.0247 / 0.3061 / 0.2545. The ⁷Li α tag is optics-blind at 0.966 / 0.974 / 0.979
(YR) against 0.981 / 0.993 / 0.994 (tagging).

## 6. HepMC3, and handing events to the T2 (PYTHIA) tier

```cpp
#include "lipolgen/hepmc_writer.hpp"
HepMC3Writer w("out.hepmc");
p.for_each([&](const Event& ev) { w.write(ev); });
w.close();
```

The T0 hadronic system is written as **documentation status 3** (pdg 92), so it
is never double-counted against real T2 hadrons; a read-back conservation check
sums status 1 plus that one status-3 pdg-92 particle.

The T2 hook has the signature of `PythiaBridge::hadronize`:

```cpp
PythiaBridge bridge(p.beam_config());
cfg.hadronizer = [&bridge](Event& ev, Rng& rng) { bridge.hadronize(ev, rng); };
```

It is called once per finished T0/T1 event with the event's own counter-based
stream, after every particle of those tiers is in place. **That is the whole
binding on every channel except coherent** — a tagged event at the default
`Tier::T1` names its `Role::StruckNucleon`, so the bridge uses it verbatim and
the whole record conserves to the numerical floor (measured 1.4 × 10⁻¹³
relative, charge exactly 0, and no "no surrogate" tail at all). The bridge's
`Role::StruckCluster` branch is deprecated: it warns once and counts into
`PythiaBridgeStats::n_cluster_fallback`, and only a `Tier::T0` run or a
hand-built record can reach it.

The cluster is still on the record as documentation, and

```cpp
StruckCluster sc;
struck_cluster_of(ev, *p.tagged_channel(), sc);
sc.p_per_nucleon_eff;   // P_X,z / A_partner
sc.virtuality;          // M_X^2 - m_free^2 < 0
```

rebuilds it from a finished event; `kin.alpha_s` / `kin.pt_s` are the
light-front variables in the ion rest frame. A hadronizer used with
`for_each(sink, nthreads)` must be re-entrant.

An **inclusive** event names its struck nucleon (`Role::StruckNucleon`, drawn
`Z F2p : N F2n`, on shell at `M_NUCLEON`), so the bridge uses it verbatim and
its own implicit-target fallback never fires on a `Pipeline` event. Driving the
bridge directly without one falls back to
`NucleonChoice::ByStructureFunctions`, which applies the same rule on
`PythiaBridgeOptions::f2_source` — hand it the kernel's own backend to keep the
two draws consistent. A **coherent** event carries `Role::Pomeron` instead of
a struck nucleon and hadronizes on the bridge's third, Pomeron-beam instance
(§4; `docs/PYTHIA_BRIDGE.md` §12).

The hard-process flavour is offered with probability `e_q² x f_q(ζ_q, Q²)` at
each flavour's OWN `ζ_q = (Q² + m_q²)/(P_A⁺ q̃⁻)`, and a flavour whose ζ_q has
run past 1 is not offered at all — `PythiaBridgeStats::n_flavour_dropped`
counts those, and they cost no `pythia.next()` retries.

## 7. Checks you get for free

```cpp
Vec4   r = momentum_residual(ev);   // zero to rounding, per the channel's balance
Vec4   s = momentum_scale(ev);      // what entered that balance
double q = charge_residual(ev);     // exactly zero
```

## 7a. What else the record carries

`Kinematics` stores two blocks a consumer would otherwise have to rebuild:

```cpp
ev.kin.cell;        // accepted-cell index of the sampler the (x, Q2) came from
                    // -- the ION-level one (inclusive), the STRUCK-CLUSTER one
                    // (tagged), the f_coh-reweighted one (coherent).  This is
                    // polligen's Mode-W `cell` column.
ev.kin.spec_pt;     // the spectator's LAB block: boost_spectator's own numbers
ev.kin.spec_theta;  // (the coherent recoil's on that channel)
ev.kin.spec_p_lab;
ev.kin.spec_r;      // rigidity ratio vs the beam; NaN for a neutral fragment
ev.kin.spec_xl;
ev.kin.spec_kx; ev.kin.spec_ky; ev.kin.spec_kz;   // rest-frame, lab-oriented
ev.kin.phi_spec;    // LAB azimuth -- NOT kin.phi, which is the DIS azimuth
```

## 7b. Radiative corrections — a band and a background, never a shift

```cpp
cfg.rc = PipelineRc::TensorBand;             // default: Off = today, bit for bit
cfg.rc_options.delta_low_x = 0.30;           // band it: run 0.19 as well
cfg.rc_options.fq_scale = 1.0;               // band it: run 0.0 and 2.0
cfg.rc_options.qe_suppression = 1.0;         // band it: run 0.0 and 0.5
cfg.rc_options.qe_tensor_scale = 0.0;        // the POLARISED QE tail: 0 = not
                                             // priced (default).  1 is a
                                             // BORROWED magnitude, not a bound
cfg.rc_options.c0_shape = C0Shape::Ho;       // band it: run C0Shape::VmcFt too
```
```bash
python -m lipolgen.cli --rc tensor-band --events 400000
python -m lipolgen.cli --rc tensor-band --rc-delta-low-x 0.19 --events 400000
python -m lipolgen.cli --rc tensor-band --rc-c0-shape vmc-ft --events 400000
```
```python
cfg = lipolgen.make_config(events=400000, rc="tensor-band",
                           rc_delta_low_x=0.30, rc_fq_scale=1.0,
                           rc_qe_suppression=1.0, rc_c0_shape="ho")
lo, hi, tail = lipolgen.export.rc_columns(pipeline.generate(0))
```

`--rc tensor-band` (`rc.hpp`, opt-in, default off) adds **three weights per
event** and changes nothing else:

| weight | what it is |
|---|---|
| `rc_tensor_lo` = 1 − δ(x)·τ | the **low** edge of a two-sided systematic band on the *tensor part alone* of the event's rate |
| `rc_tensor_hi` = 1 + δ(x)·τ | the **high** edge of the same band |
| `rc_tail` = 1 + σ_tail/σ_Born | the ⁶Li **radiative-tail background**: the elastic tail (with its small tensor part) plus the *unpolarised* quasi-elastic tail |

with `τ = W_tensor/W` the rank-2 fraction of the event's own density —
`(P_zz/2)A_zz / [1 + (P_zz/2)A_zz]` at θ_S = 0 — and δ(x) log-linear between
two anchors, `δ = 0.30` at `x = 0.01` and `δ = 0.015` at `x = 0.16`, clamped
outside them.

**They are NOT on `Event::weight`.** FSI multiplies the nominal weight
because FSI is a *correction to the model*; these are not. The band is a
*systematic variation* and the tail is a *background*, and both must leave
the Born sample alone. They travel in `Event::rc_weights` under their own
HepMC3 names and their own npz columns, and an analysis multiplies one in on
purpose:

```python
w_hi = cols["weight"] * cols["rc_tensor_hi"]     # the +delta edge
w_bg = cols["weight"] * cols["rc_tail"]          # Born + radiative tails
```

`rc_tensor_lo + rc_tensor_hi == 2.0` bit-for-bit on every event: they are the
two **signs** of one δ, not two different band edges.

**No four-vector moves and no random number is consumed**, so an
`--rc tensor-band` run is bit-for-bit an `--rc off` run in every particle,
every kinematic label, `Event::weight` and the RNG stream — proved event by
event, and file by file, in `tests/test_rc_pipeline.cpp` (T6) and
`python/tests/test_rc.py`. The structural guarantee is that
`RcModel::fill(Event&) const` takes no `Rng&`. With `--rc off` the npz key
set and the HepMC3 weight names are **exactly** today's, so an `--rc off`
file is byte-identical to one written before `rc.hpp` existed.

**Per channel** (`RcModel::applies()` decides; the run prints why, and
nothing is ever refused so a channel scan need not special-case `--rc`):

| channel | band | tail | why |
|---|---|---|---|
| inclusive | full, at every θ_S | full, at every θ_S | the home case |
| every `Tagged*` | full (τ from the cluster density), **clamped** | **≡ 1**, and that is HALF a kinematic fact | the **elastic** recoil sits at x_L = 1, inside the 10σ beam envelope, so the tag itself vetoes it — a fact. The **quasi-elastic** tail is **not** vetoed and its absence is an **omission**: see the box below |
| coherent ⁶Li | **≡ 1** | **≡ 1** | its tensor dependence is entirely azimuthal, and **no RC treatment exists for a φ-dependent tensor observable**; the elastic point (M_X = 0) is already outside the channel |

**Two things the "full" in that table does not say.**

* **On ⁷Li the CLI reaches the band only at `--pe 0`.** `--rc tensor-band`
  needs an **unpolarised beam** (`RcModel` throws on any category with
  `λ_e·P_e ≠ 0` — the whole `A_zz` programme assumes one). The one CLI route
  is therefore

  ```bash
  lipolgen-run --isotope 7Li --channel tagged-7Li-alpha \
               --plan helicity-flip --pe 0 --rc tensor-band --events 300
  ```

  where `λ_e·P_e = 0` satisfies the check and the run writes
  `rc_tensor_lo`/`rc_tensor_hi`/`rc_tail`. **Every other CLI plan is
  refused**: the helicity plans (`apar`, `helicity-flip`) at `P_e ≠ 0` by
  `RcModel`, and the spin-1 tensor plans (`azz`, `tensor-thirds`,
  `transverse-tensor`, `cos2phi`, `tensor-flip`) by `Pipeline` with "run-plan
  spin 1.0 != channel ion spin". All eight were checked. An explicit
  (P_z, T) `J = 3/2` fill — a *polarised* ⁷Li target with an unpolarised beam,
  which is what the physics case actually wants — still needs the API, because
  no CLI plan builds one:

  ```python
  cats = [_l.SpinCategory("m32", 1.5, [0.5, 0.0, 0.0, 0.5]),
          _l.SpinCategory("m12", 1.5, [0.0, 0.5, 0.5, 0.0])]
  plan = _l.RunPlan(cats, 0.0, 0.0, 0.6)          # unpolarised beam, J = 3/2
  p = lipolgen.Pipeline(make_config(isotope="7Li", channel="tagged-alpha",
                                    rc="tensor-band"), plan)
  ```

  (`python/tests/test_rc.py` runs exactly this, asserts that
  `helicity-flip` at `pe = 0` **does** build a working band, and re-checks
  that every CLI plan at `P_e ≠ 0` is refused — so the claim cannot rot.)

> **`rc_tail ≡ 1` on a tagged channel is half a kinematic fact — read this
> before treating a tagged run as tail-free.** The ELASTIC half is a fact: the
> intact ion recoils at `x_L = 1`, inside the 10σ beam-exclusion envelope,
> while the tag looks at `x_L ≈ A_spec/A_beam`. The QUASI-ELASTIC half is an
> **omission**, and the veto argument does not cover it:
>
> * quasi-elastic knockout removes **one** nucleon; the A−1 remnant of either
>   lithium channel is **unbound** — ⁵Li and ⁵He are resonances above the
>   α + N threshold with no particle-stable state at all — so it breaks up and
>   its α emerges at
>   `x_L ≈ (4/5)(5/6) = 2/3` — **the tag window itself**;
> * in the α + d picture, if the struck nucleon is one of the embedded
>   deuteron's two the α is a **true spectator** carrying the same `n_M(k, c)`
>   the tagged Born is built on. On the deuteron control the quasi-elastic tail
>   simply **is** elastic `e`–`n` scattering with a spectator proton, the
>   classic spectator-tagging background;
> * so the tag does not suppress it in the **ratio** `rc_tail` is: the same
>   spectator density and Roman-Pot acceptance multiply the tail and the tagged
>   Born, and both select the **2 of 6** nucleons inside the deuteron (a
>   nucleon knocked out of the **α** destroys the α and *is* vetoed).
>
> The omitted dilution is therefore of the **same order as the inclusive
> quasi-elastic one** — 22 % / 73 % / 99.9 % of the inclusive tail at
> x = 0.01 / 0.10 / 0.30 — not a negligible one. That is an
> **order-of-magnitude argument, not a computed number**; pricing it needs a
> **tagged** Born denominator and the tag acceptance folded into Eq. (44), plus
> the **cluster-elastic** `e + A → e' + γ + d + α` that Eq. (44)'s free-nucleon
> sum does not contain at all. `RcModel::exclusion_reason()` says all of this
> in the run banner and in `meta["rc_exclusion_reason"]`, so a tagged
> `rc_tail == 1` cannot be read as a veto on the whole tail.
> (`docs/open_items/run_2026-09-03/phase_B_numbers.md` §B4.)

* **On the tagged channels the band is CLAMPED, and it has to be.**
  `τ_tag = 1 − n̄(k,c)/n_M(k,c)` is unbounded: wherever the event's own
  `n_M` is near a node of the M-dependent spectator density — the ⁶Li `M = 0`
  density has them — the ratio blows up (measured **|τ| up to 30.7**).
  Unclamped, a 20 k-event tagged-alpha run published `rc_tensor_hi` down to
  **−1.79** (⁷Li: **−8.35**), i.e. **negative weights in the npz**.
  `RcOptions::band_tau_max` (default 1.0) holds every published edge inside
  `[1 − δ, 1 + δ]`; the clipped fraction is **1.2 % of ⁶Li tagged-alpha
  events and 2.6 % of ⁷Li** (re-measured 2026-09-15; the ⁶Li 0.6 % predates
  the S–D fix `a7b3d18`), is printed by the run and lands in
  `meta["rc_clipped_band_event_fraction"]` and per event in
  `Event.rc_clipped`. `n_M → 0` is exactly where "rescale the tensor part by
  δ" stops meaning anything, so the clamp is a **choice**, not a fix.

**The mandatory band runs — the 20–40 mb rule of `--fsi-sigma-mb`, applied
here.** Never quote one row alone. Run each knob at both/all its ends and
quote the lo/hi envelope on `A_zz`:

| knob | rows to run | what it prices |
|---|---|---|
| `--rc-delta-low-x` | **0.19 and 0.30** (and, since 2026-09-06, **0.266 and 0.113** as priced alternatives) | 0.30 is the conservative end of Gakh–Shekhovtsova's 10–30 %; 0.19 is the residual HERMES actually achieved at its lowest-x bin. **Which way the default errs:** the 0.30 is that paper's panel value **carried upward in x**, and **0.113 is what the panel reads at x = 0.00966, the x nearest the 0.01 anchor** (0.266 at its own bottom, x = 0.00226) — the shipped anchor is **×2.65** that reading. Priced, not adopted (registry row 17, `open_items/run_2026-09-06/phase_A_numbers.md` §A2): band half-width on A_zz at Q² = 5, P_zz = +1, ⁶Li config 1 = **1.358472e−04 / 1.204512e−04 / 5.116913e−05** at x = 0.01, **1.459067e−04 / 1.308566e−04 / 6.313127e−05** at x = 0.063 (where the band peaks on all three), **9.680016e−05 / 8.798804e−05 / 4.833349e−05** at x = 0.10 and **1.920691e−05** at x = 0.16 on all three — the high anchor pins x ≥ 0.16, so the whole choice lives below it. It moves **two columns and two of the 60 `meta` keys** (`rc_delta_low_x` and the `knob_provenance` row that records it) and **not** the tagged band's clipped-event fraction (245 of 20 000 = 1.2250 % on all three, seed 1 — re-measured 2026-09-15; 124 = 0.62 % predates the S–D fix `a7b3d18`: `clamp_tau` clips \|τ\| and δ never enters it) |
| `--rc-a-transfer-frac` | **0 (default), 0.5, 1** | the **A = 2 → A = 6 TRANSFER** of the band. Every number δ(x) interpolates between is a **deuteron** number — HERMES's measured low-x residual, Gakh–Shekhovtsova's 10–30 % (itself a Q² = 0.1 GeV² figure — the only panel of its Fig. 2 inside the x ∼ 10⁻³–10⁻² the sentence quotes it for), E12-13-011's 1.5 % — and **no A > 2 tensor RC calculation exists at all**, so the default band silently assumes the deuteron *fractional* RC transfers to ⁶Li exactly. This adds `f·δ(x)` in quadrature: `δ_eff = δ√(1+f²)`, i.e. the band widens by **1× / 1.118× / 1.414×**. Half-widths on A_zz at x = 0.01, Q² = 5: **1.358e−04 / 1.519e−04 / 1.921e−04**. It moves the **band only** — not `rc_tail`, not τ — and **no measurement prefers any value**: it is a price tag, not a correction. `meta["rc_a_transfer_frac"]` records it |
| `--rc-fq-scale` | **0, 1, 2** | ±100 % on the ⁶Li quadrupole form factor. σ^el_T is **quadratic** in it, so this band must be **RUN, never rescaled** from one row — the two edges are not symmetric about the nominal, and at `x = 0.01, Q² = 5` they even bracket a **sign change** of ΔA_zz |
| `--rc-tail-tensor-scale` | **0.5, 1, 2** | the η·F_m² tensor sector, which `--rc-fq-scale` does **not** span |
| `--rc-c0-shape` | **`ho` and `vmc-ft`** | the ⁶Li **C0 (monopole) SHAPE**, shared by F_c and F_q. `ho` is the unfitted harmonic oscillator every published number was made with; `vmc-ft` is the j₀ transform of the committed ANL VMC point-proton density, r-rescaled so that **⟨r²⟩_point, F_c(0) = 3 and F_q(0) = −65.914 are identical on the two edges** (T11 gates all three on both). It is a *shape*, not a multiplier: it cannot be rescaled out of one run, and no two-parameter oscillator can be refitted to reach the other edge while holding ⟨r²⟩. It **flips the sign of the tensor fraction of the elastic tail at x = 0.1** — see the box below. Cost: `RcModel` construction 0.11 s → 1.29 s |
| `--rc-tail-model` | **`t-peak`, `t-peak+ll` and `polrad-full`** | **WHICH PEAKS OF THE TAIL ARE IN IT.** `polrad-full` (2026-09-06) is POLRAD **Eq. (18) + Appendix B + Eq. (A.4)** — ONE exact τ_A quadrature with all three peaks **and** the s-/p-peaks' own Eq. (A.4) tensor content, so `--rc-sp-tensor-scale` is refused on it (the term RAN). **It is NOT inside the t-peak pair**: between the two edges on **1725 of 3051** accepted cells (56.5 %), covering a median **0.6555** of the gap **over the 3027 cells whose gap is nonzero** (0.6620 over all 3051; the other 24 have `t-peak+ll` = `t-peak` exactly, so no fraction), with the **ratio of the σ-weighted mean shifts** at **0.428** — a ratio of means, **not** a σ-weighted (still less an event-weighted) mean of the per-cell fractions, which is **6.483** — and **above both** in the Q² ≥ 20, y ≤ 0.9 window. Whole-run mean `rc_tail` (200 k, seed 1234) **at the default fill P_z = 0.7**: **1.021778527 / 1.040274188 / 1.029702912**; at **P_z = 0** (the plan both test suites use) **1.021836305 / 1.040368933 / 1.029775347** on the same build. Still **not** checked against [MT69] or any external exact tail — only POLRAD-internally (`Eq. (18)/Eq. (38) → 1` as `x_A → 0`, 1.00230 unpolarised / 1.00313 tensor at x_A = 0.003) and against the leading log. Costs ~10 s of table build and needs `n_eta ≥ 64`. `t-peak` is the default and the **lower** edge: POLRAD Eqs. (37)–(39), (43), one peak of a three-peak object. `t-peak+ll` adds the leading-log s- and p-peaks and is the **upper** edge. Their sum is a **stated model of mixed approximation orders** (an η_A quadrature plus a single-z leading log), good to ~5–10 %, with an uncancelled soft `1/(1−z)` as `y → 0`; and it has **no tensor s/p partner**, so it **lowers the tensor fraction** of the tail (`r_T/r_U` ×0.66139 / ×0.0031829 / ×6.6076e−05 at x = 0.01 / 0.10 / 0.30, Q² = 5) — the tensor part of those peaks is, since 2026-09-06, **bounded, not computed** by `--rc-sp-tensor-scale`, and that bound is **empty at those three points** (recovers 0.0 % of the collapse, bit-identical) and at most **21.72 %** anywhere on the grid, with the quasi-elastic s/p column bounded by **nothing** — see §7b. Whole-run mean `rc_tail` on 200 k inclusive ⁶Li config-1 events, **at the default fill P_z = 0.7**: **1.021778527 → 1.040274188**; at **P_z = 0**, **1.021836305 → 1.040368933** (*corrected 2026-09-15: the 2026-09-06 note that the 1.021836 → 1.040369 pair "is pre-Phase-A and does not reproduce" is **withdrawn** — it is this build's own P_z = 0 answer, re-measured*). In the `Q² ≥ 20 GeV²`, `y ≤ 0.9` window the two edges agree to **+0.61 % event-weighted at P_z = 0, +0.62 % at P_z = 0.7** but **not cell by cell**: 331 of 1356 accepted cells (24.4 %) differ by more than 1 %, worst **×6444** at `x = 0.79`, `y = 0.0088`, and they differ by 59 % at `y = 0.985`. Per-cell agreement is ≤ 0.55 % only for `0.15 ≤ y ≤ 0.7` |
| `--rc-qe-suppression` | **0, 0.5, 1** | a flat multiplier on the quasi-elastic tail, **on top of** the Pauli suppression below. `rc_tail` is quasi-elastic-**dominated** at every `x ≳ 0.03` (73 % of the tail at `x = 0.1`, **99.9 %** at `x = 0.30`), so 0 is never a small variation |
| `--rc-qe-tensor-scale` | **0 (default) and 1 — and read the caveat before either** | the **POLARISED** quasi-elastic tail, which is otherwise treated as exactly tensor-blind on the piece that is 73 % of `rc_tail` at x = 0.1 and **99.9 %** at x = 0.30. Nobody has computed it: POLRAD has no tensor partner to Eq. (44) and no such calculation exists for an A = 6 spin-1 nucleus. At **1** the quasi-elastic tail is lent the **elastic** tail's own σ^el_T/σ^el_U — a **BORROWED MAGNITUDE, not a derived bound** (a *coherent* nuclear quadrupole fraction on an *incoherent* nucleon process). ⁶Li's elastic tensor fraction is anomalously small for a reason the quasi-elastic piece has no reason to share, so **1 may be ~10² too small at x ≤ 0.1**, and the SIGN it inherits is meaningless — read the **magnitude**. It is **exactly linear**, so one run rescales to any value. Measured at Q² = 5: ΔA_zz moves by −1.162e−07 / +3.018e−10 / +1.796e−09 at x = 0.01 / 0.10 / 0.30, i.e. **0.086 % / 0.0003 % / 0.25 %** of the band half-width there |
| `--rc-sp-tensor-scale` | **0 (default) and 1 — and read where the bound is EMPTY before quoting either** | the **TENSOR fraction of the leading-log s-/p-peaks**, which `--rc-tail-model t-peak+ll` otherwise adds to the **unpolarised** numerator alone, so that edge lowers `r_T/r_U` by **×0.66139 / ×0.0031829 / ×6.6076e−05** at x = 0.01 / 0.10 / 0.30, Q² = 5 purely by growing the denominator (`r_U` ×1.512 / ×314.18 / ×15134). **Refused unless `--rc-tail-model t-peak+ll`** — the t-peak-only tail computes no s/p peaks, so a price there would be recorded without a single operation behind it. At **1** the **elastic** s-/p-peaks are lent the elastic t-peak's own σ^el_T/σ^el_U. **IT IS A BOUND WITH NO DERIVATION**, and it borrows *less* than `--rc-qe-tensor-scale` (the same coherent ⁶Li vertex, a different photon topology) — but **it is EMPTY at the three standard points**: Q′²_s = 4.3728 / 4.9382 / 4.9801 GeV² there, where F_c = −3.75e−45 / −6.44e−51 / −2.41e−51 against **+2.99260** at the t-peak's own t_min = 8.7304e−05, so `u_sp/σ^el_U` = 5.25e−79 / 2.07e−86 / 4.06e−83 and **ΔA_zz moves by exactly 0 at scale 0, 0.5 and 1**, recovering **0.0 %** of the collapse. It bites on **77 of 3051 accepted cells** (x ≤ 7.94e−03, y ≥ 0.366), reaching **582.9 %** of the band half-width at x = 4.169e−04, y = 0.9692 (`--rc-qe-tensor-scale 1` gives 328.4 % there) and recovering at most **21.72 %** of the collapse. Exactly linear, so one run rescales. The **quasi-elastic** s/p column is bounded by **neither** scale |
| `RcOptions::qe_kf_gev` (API) | **0.169 (default), 0.221, 0** | POLRAD Eq. (44)'s `S_E`/`S_M`, the de Forest–Walecka Fermi-gas factor `S(q) = (3/4)(q/k_F) − (q/k_F)³/16` below `q = 2k_F`, exactly as POLRAD's `ffquas` codes it. **On by default** at ⁶Li's measured `k_F` (Moniz *et al.*, PRL **26** (1971) 445). It cuts the QRT to **0.47** at `x = 0.01` and **0.87** at `x = 0.1`; `0` is the unsuppressed edge v0 shipped |

> **The C0 shape band decides a sign, so read this before quoting `σ^el_T`.**
> `(1/6)σ^el_T/σ^el_U` at Q² = 5 GeV², `ho` / `vmc-ft`:
>
> | x | 0.01 | 0.03 | 0.10 | 0.30 |
> |---|---|---|---|---|
> | `ho` | −5.0943e−04 | −7.3583e−04 | **+1.5595e−04** | +1.3287e−02 |
> | `vmc-ft` | −5.4705e−04 | −8.0426e−04 | **−4.0702e−05** | +1.5556e−02 |
>
> The **sign change with x survives** the band — both edges are negative at low
> x and positive at 0.30 — but **where** it happens does not, and **at x = 0.10
> the sign is indeterminate**: publish it as such, never as +1.6e−04. Adding
> `--rc-fq-scale 0 … 2` widens x = 0.10 to −9.0e−05 … +4.0e−04, still spanning
> zero on both edges. `σ^el_U` rises by ×1.011 / ×1.022 / ×1.104 / ×3.03 over
> the same x, the quasi-elastic tail is **bit-identical** between the edges, and
> `ΔA_zz` moves by −8 % / +0.6 % / +37 % at x = 0.01 / 0.10 / 0.30 — three
> orders below the band, so the *ordering* of the systematics is unchanged.
> `meta["rc_c0_shape"]` records which edge ran. Everything measured, with the
> high-q continuation and the two policy overrides argued:
> `docs/open_items/run_2026-09-03/phase_B_numbers.md` §B1.

At `x = 0.01, Q² = 5 GeV²` (config 1, `--plan tensor-thirds --pzz 0.6`) the
band is by far the largest of these: half-width **`1.358472e−04`** on `A_zz` at
δ_low = 0.30, against `8.0e−07` from the whole `fq_scale`
band, `2.9e−08` from the whole C0 shape band and `1.769797e−07` from the entire
tail — a factor **768** over the tail. **The band — the unapplied
lepton-vertex correction — is what `--rc tensor-band` is for; the tail is a
bookkeeping item at EIC energies and the headline at a fixed target.** Every
number is in `docs/OPEN_ITEMS_SOLUTIONS.md` §9.

> **Corrected here 2026-09-06, and one entry deliberately left flagged.**
> This paragraph read `4.4e−04` (and `2.8e−04` at δ_low = 0.19) against
> `3.5e−07` from the tail. Both carried the **×3.253983 deuteron-b₁ `A_zz`**
> withdrawn on 2026-09-04
> (`open_items/run_2026-09-03/phase_B_numbers.md` §B3.2); with this
> generator's own ⁶Li b₁ the half-width is **`1.358472e−04`** (and
> **`8.603659e−05`** at δ_low = 0.19, re-measured here) and the whole tail
> **`−1.769797e−07`**, so the factor is **768**, not 1300.
> **The other two entries — `8.0e−07` from `fq_scale` and `2.9e−08` from the
> C0 shape — are `OPEN_ITEMS_SOLUTIONS.md` §9's own un-re-measured rows and
> still carry the withdrawn `A_zz`**; they are quoted unchanged and are not
> to be read as corrected. The **ordering**, which is all this paragraph
> exists to state, does not depend on them: both sit three to four orders
> below the band on any reading. Registry row 17's own price — three anchors,
> four x — is `open_items/run_2026-09-06/phase_A_numbers.md` §A2.

**Honest flags, to repeat wherever any of this is quoted:**

* `δ_low = 0.30` is the **size of a correction this generator does not
  apply**, taken as a 1σ band — *not* a measured residual. Its source
  (Gakh–Shekhovtsova, [hep-ph/0403262](https://arxiv.org/abs/hep-ph/0403262))
  has **zero INSPIRE citations**. And it is that source's panel value
  **carried upward in x**: the panel spans x = 0.00226–0.00966 at Q² = 0.1,
  and **0.113 is the value it actually reads at x = 0.00966 — the x nearest
  the 0.01 anchor** (0.266 at its own bottom). The shipped anchor is **×2.65**
  that reading, i.e. it errs **wide**. Both alternatives are **priced, not
  adopted** (registry row 17,
  `open_items/run_2026-09-06/phase_A_numbers.md` §A2).
* Everything cited is **deuteron**. Whether the deuteron's *fractional* RC
  transfers to ⁶Li is untested, and is the largest unquantified assumption
  in the band.
* The ⁶Li form-factor shape parameters are **unfitted starting values**, not
  a fit (the elastic data are not in this repository in machine-readable
  form). The *normalisations* are the **measured** moments — `μ = +0.822047
  μ_N`, `Q = −0.0818 fm²` — and never VMC, whose `Q(⁶Li) = −0.23(9) fm²` is
  3× the measured one Since 2026-09-04 the **C0 shape** is banded rather than
  left to look settled (`--rc-c0-shape`, above); a band is a **price tag, not a
  fit**, and Q1 stays open. Its `vmc-ft` edge does read the VMC density — for
  the **monopole shape only**, never for a tensor normalisation, which is what
  `OPEN_ITEMS_SOLUTIONS.md` §5 forbids; and it makes
  `data/vmc/density/li6.density` an *input*, where `design_G_cluster_config.md`
  §3.2 had declared it the independent validation target of the cluster model's
  T6. Both overrides are written down in `include/lipolgen/rc.hpp` beside the
  constants and in `phase_B_numbers.md` §B1.7 — `src/core/cluster_config.cpp`
  still opens the file nowhere, so T6 itself stays independent.
* The ⁶Li form-factor **dip location is a model number on either edge**:
  q₀ = 3.0999 fm⁻¹ ⇒ `|t| = 0.3742 GeV²`, gated over [2.9, 3.3] fm⁻¹ by T11 —
  a starting guess no data in this repository can refit — while the committed
  VMC point-proton density has no C0 zero below q ≈ 4.3 fm⁻¹ at all, and the
  `vmc-ft` edge accordingly has none anywhere. The uncited "`|t| ≈ 0.31 GeV²`"
  that `PHYSICS_CHANNELS.md` carried until 2026-09-04 is withdrawn.
* **No RC calculation exists for a tagged tensor asymmetry**, so the tagged
  band is a defensible but **uncited extrapolation**; and none exists for any
  φ-dependent tensor observable at any axis, which is why the `Δ` cos 2φ
  sector and the coherent channel carry no band at all.
* The **polarised** quasi-elastic tail is **not computed anywhere** — POLRAD
  has no tensor partner to Eq. (44) and no such calculation exists for an
  A = 6 spin-1 nucleus — and it is **neglected by default**, citing Z.-L. Zhou
  *et al.*, PRL **82** (1999) 687 (a *deuteron* statement, and a paper not in
  this repository). Since the quasi-elastic term is **22 % / 73 % / 99.9 %**
  of `rc_tail` at x = 0.01 / 0.10 / 0.30, that zero is the largest unpriced
  piece of the tail, and `--rc-qe-tensor-scale` (default **0.0**, the shipped
  tensor-blind tail bit for bit) prices it by lending the quasi-elastic tail
  the **elastic** tail's own σ^el_T/σ^el_U. **That is a borrowed magnitude,
  not a derived bound**: it puts a coherent nuclear quadrupole fraction on an
  incoherent nucleon process, ⁶Li's elastic tensor fraction is anomalously
  small for a reason the quasi-elastic piece does not share, and the SIGN it
  inherits is meaningless. Read `RcOptions::qe_tensor_scale` in `rc.hpp`
  before quoting a number from it; measured sizes in
  `docs/open_items/run_2026-09-03/phase_B_numbers.md` §B3.
* **`rc_tail` has THREE tail models, the two t-peak ones are a PRICE RANGE
  and NOT a confidence interval, and none of the three is "the" radiative
  tail.** `--rc-tail-model t-peak` (the **default**, bit for bit every
  published number) is ONE PEAK of the elastic and quasi-elastic tails and its
  absolute normalisation is not validated against any exact tail: it is a
  **lower bound on the dilution**. `--rc-tail-model t-peak+ll` adds the
  leading-log s- and p-peaks of the same two unpolarised observables
  (`ll_peaks_spin1` / `ll_peaks_qe`, promoted out of the test binary in the
  2026-09-03 run) and is the **upper edge**. `--rc-tail-model polrad-full`
  (2026-09-06) is POLRAD **Eq. (18) + Appendix B + Eq. (A.4)**: ONE exact τ_A
  quadrature carrying all three peaks, with the s-/p-peaks' own Eq. (A.4)
  **tensor** content. Run all three; quote all three.

  **AND THE EXACT ONE IS NOT INSIDE THE OTHER TWO.** Measured over the
  sampler's own 3051 accepted ⁶Li config-1 cells: `polrad-full` lies between
  the two t-peak edges on **1725 (56.5 %)** and **outside on 1326 (43.5 %)**;
  it covers a **median 0.6555** of the `t-peak → t-peak+ll` gap **over the
  3027 cells whose gap is nonzero** (0.6620 over all 3051 — on the other 24
  `t-peak+ll` equals `t-peak` exactly and the fraction is undefined), and the
  **ratio of the σ-weighted mean shifts** is **0.428** (0.427958 unclipped,
  0.427368 at the shipped `tail_max = 10`) — a **ratio of means**, *not* a
  σ-weighted mean of the per-cell fractions, which is **6.483**; and in the
  `Q² ≥ 20 GeV²`, `y ≤ 0.9` window it is **above
  both** (8.815031e−03 against 8.719649e−03 and 8.773752e−03, at the default
  fill P_z = 0.7). Whole-run mean
  `rc_tail`, 200 k events, seed 1234, **at P_z = 0.7**: **1.021778527 →
  1.040274188 → 1.029702912**; at **P_z = 0** (`tensor_thirds_plan(0.0, 0.6)`,
  the plan both test suites use) **1.021836305 → 1.040368933 → 1.029775347**
  on the same build. The cell-σ-weighted census above is P_z-free; an
  event-weighted mean is not, because the fill plan decides which events the
  sampler draws (re-measured at both, 2026-09-15).

  **What `polrad-full` is still NOT checked against: Mo–Tsai, or any external
  exact tail.** `[MT69]` is not in this tree and no number from it is quoted
  anywhere. What is checked is **POLRAD-internal** — with a form factor dead
  at the s-/p-peak vertex, so only the t-peak survives, `Eq. (18)/Eq. (38)` is
  **1.00230** (unpolarised) and **1.00313** (tensor) — the F_m-only sector's ratios; the spin-0 unpolarised ratio at the same x_A is 1.00492 and the F_q tensor ratio 1.01872 (§B2.3) at `x_A = 0.003` and tends
  to 1 as `x_A → 0` — the `Q_N = 0` Rosenbluth limit of Eq. (A.4) (T9,
  un-skipped for this), and the leading-log fallback. It costs **~10 s** of
  tail-table build and needs `rc_options.n_eta ≥ 64`; both a smaller `n_eta`
  and a platform whose `long double` is no wider than `double` are **refused**,
  never silently degraded. Details, and the FIVE transcription defects found in
  `polrad2t.tex`'s Appendix B on the way, in
  `docs/open_items/run_2026-09-06/phase_B_numbers.md` §B2.

  **The upper edge is a STATED MODEL, not a controlled expansion.** It sums
  POLRAD's η_A quadrature and a single-z collinear leading log, so it is
  accurate to the worse of the two (**~5–10 %**, the size of `1/ln(Q²/m_e²)`),
  and its radiator carries an **uncancelled soft `1/(1−z)`** that overshoots
  as `y → 0` — where `z_s → 1` the radiator exceeds 1 and single emission is
  no longer the right expansion. The absolute contribution there is small
  (`< 4e−4` of the Born at `x = 0.74`, `y = 0.007`) but the model is not
  trustworthy in that corner.

  **NEITHER t-PEAK MODEL carries a tensor s/p peak**, because **Eq. (38)**
  supplies none and a leading log cannot be given one without a second
  definition of a physics number (`docs/CONVENTIONS.md`). **Eq. (18) DOES**,
  and `polrad-full` computes it — the sentence *"POLRAD supplies no tensor s/p
  peak"*, which stood here until 2026-09-06, was true of Eq. (38) and false of
  the paper. On the two t-peak models the s+p therefore enter the
  **unpolarised** numerator only, so `t-peak+ll` **lowers the tensor fraction
  of the tail** —
  measured at Q² = 5 GeV², ⁶Li config 1, production grid, 2026-09-06: `r_T/r_U`
  falls by **×0.66139 / ×0.0031829 / ×6.6076e−05** at x = 0.01 / 0.10 / 0.30,
  because `r_U` grows ×1.512 / ×314.18 / ×15134 while the tensor numerator does
  not move at all.

  Since 2026-09-06 the tensor part of those peaks is **bounded, not computed** —
  and the number belongs with the word. `--rc-sp-tensor-scale`
  (`RcOptions::sp_tensor_scale`, default **0.0**, **refused** unless
  `--rc-tail-model t-peak+ll`, exactly linear so one run rescales) lends the
  **elastic** s-/p-peaks the elastic t-peak's own σ^el_T/σ^el_U: scale 1 is the
  sentence *"the s/p tensor fraction equals the elastic t-peak's"*. **It is a
  bound with no derivation.** It borrows *less* than `--rc-qe-tensor-scale` —
  the same coherent ⁶Li vertex reached by a different photon topology, not a
  coherent ratio lent to an incoherent process — but it is **EMPTY at the three
  standard points**: the s/p elastic vertex sits at Q′²_s = **4.3728 / 4.9382 /
  4.9801 GeV²**, where F_c = −3.75e−45 / −6.44e−51 / −2.41e−51 against
  **+2.99260** at the t-peak's own t_min = 8.7304e−05 GeV², so `u_sp/σ^el_U` =
  5.25e−79 / 2.07e−86 / 4.06e−83 and scale 1 is **bit-identical to 0** there.
  **And now that there is a computed answer to score it against, the bound was
  not even one-sided**: at the same three points `polrad-full` gives
  `r_T/r_U` **×0.79395 / ×0.0020032 / ×0.00022591** of the t-peak against
  `t-peak+ll`'s ×0.66139 / ×0.0031829 / ×6.6076e−05 — so at x = 0.01 and 0.30
  the bound pointed the right way and stopped far short, and at x = 0.10 it
  pointed the **wrong** way. `--rc-sp-tensor-scale` is **refused on
  `polrad-full`**, for the opposite reason it is refused on `t-peak`: the term
  it stands in for **ran**.
  ΔA_zz moves by **exactly 0** at x = 0.01 / 0.10 / 0.30, Q² = 5 for scale 0,
  0.5 **and** 1, against `--rc-qe-tensor-scale 1`'s **0.0855 % / 0.00031 % /
  0.246 %** of the band half-width on the same edge. It bites on **77 of 3051
  accepted cells** (x ≤ 7.94e−03, y ≥ 0.366), where at scale 1 it reaches
  **582.9 %** of the band half-width (x = 4.169e−04, Q² = 1.608, y = 0.9692,
  against `--rc-qe-tensor-scale 1`'s 328.4 % there) — and that is the corner
  where `t-peak+ll` is itself least trustworthy, so a large price there is not a
  licence to quote it. The distribution over those 77 is extreme: **3** cells
  exceed 100 % of the band, 6 exceed 10 %, 8 exceed 1 %, 18 exceed 0.1 %, and
  the median is **2.8e−05 %**. Whole-run mean `rc_tail` (200 k, seed 1234, config 1)
  moves **1.040274188 → 1.040274193** at scale 1 and → 1.040274711 at scale 100 (event-weighted at the CLI default P_z = 0.7).

  **What is still bounded by nothing:** the **QUASI-ELASTIC s/p column**
  (`TailTriple::qe_sp`). `--rc-qe-tensor-scale` multiplies σ^q_U alone and
  `--rc-sp-tensor-scale` the elastic s/p column alone, so the column that
  *survives* where the coherent one dies keeps a tensor part of exactly zero —
  and it is essentially the whole of the collapse at x ≥ 0.10. It is now the
  largest exactly-zero tensor term in `rc_tail`
  (`open_items/run_2026-09-06/phase_B_numbers.md` §B1).

  **Where POLRAD §2.1.3 B's "the s- and p-peaks are suppressed" holds, and
  where it does not** (`tests/test_rc.cpp` T8(c)/T8(d)/T8(d)(i),
  `phase_B_numbers.md` §B2). It is an **event-weighted** statement about this
  generator's bulk and it is **false cell by cell**; the two must not be
  quoted for each other.

  * **Event-weighted — it holds, and SAY WHICH P_z.** For ⁶Li at EIC config 1
    restricted to `Q² ≥ 20 GeV²` and `y ≤ 0.9`, at the CLI's **default fill
    P_z = 0.7** (`tensor_thirds_plan(0.7, 0.6)`; 5194 of 200 000 events, seed
    1234) the mean
    dilution `⟨w_tail − 1⟩` moves **8.719649e−03 → 8.773752e−03**,
    **+0.62 %**, when the s-/p-peaks are added, and **8.815031e−03**
    (**+1.09 %**) on the exact `polrad-full` tail — which is **above both
    edges**, so even this event-weighted statement does not BRACKET the exact
    answer. At **P_z = 0** (`tensor_thirds_plan(0.0, 0.6)`, the plan both test
    suites use and the plan the 2026-09-03 run published) the same build gives
    **5182** events and **8.489138e−03 → 8.540716e−03 → 8.581236e−03**
    (**+0.61 %**, **+1.08 %**). A *rate* analysis in that window is unaffected
    at the 1 % level on either.
    *(**CORRECTION, 2026-09-15 — re-measured at both P_z.** From 2026-09-06
    this bullet said the "5182 of 200 000 events … 8.48914e−03 → 8.54072e−03,
    +0.61 %" absolutes "no longer reproduce" and attributed the move to that
    run's Phase A changing the ⁶Li kernel normalisation and with it the
    sampler's cell weights. **Both claims are withdrawn.** Those digits are
    this build's own P_z = 0 answer, reproduced here to every printed digit;
    the 2026-09-06 numbers are the same run at the CLI default P_z = 0.7. The
    re-measurement had switched fill plans without saying so —
    `open_items/run_2026-09-03/phase_B_numbers.md:836-839` already tabulated
    both rows. Phase A moved nothing here. The cell-weighted
    8.31257e−03 → 8.35972e−03 pair remains withdrawn: it was never
    re-measured.)*
  * **Per cell — it fails.** Of the sampler's 3051 accepted cells, 1356 sit at
    `Q² ≥ 20` and `y ≤ 0.9`, and **331 of those (24.4 %, 28.2 % of the
    window's cross section) disagree by more than 1 %**, worst **×6444** at
    `x = 0.7943`, `y = 0.0088`, `Q² = 27.8` — where the t-peak is 0.016 % of
    the leading-log total. **318 of the 331 sit at `y < 0.1`**, i.e. the
    failure is at `y → 0`, not `y → 1`, and the `y ≤ 0.9` qualifier does
    nothing about it.
  * **The narrow per-cell claim that survives.** At `Q² ≥ 20 GeV²` **and**
    `0.15 ≤ y ≤ 0.7` no accepted cell disagrees by more than **0.55 %** (660
    cells); to `y ≤ 0.8`, **0.77 %** (725); to `y ≤ 0.9`, **2.0 %** (781).
    **Below `y = 0.15` there is no agreement statement at all.** The
    "0.16 % at `y ≤ 0.7`" this bullet used to carry was read off four table
    rows — three at `y = 0.5`, one at `y = 0.7` — and understates even the
    slice it was meant to cover.

  **Two bad corners, two different mechanisms.** At `y → 1` the s-peak beats
  `Y₊` because `z_s = (1−y)/(1−x_A y) → 0` puts the elastic vertex at
  `Q′² → 0` where the form factor is 1 (**59 %** at `x = 0.01`, `y = 0.985`).
  At `y → 0` it is the **uncancelled soft radiator**: `1 − z_s =
  y(1−x_A)/(1−x_A y)` and `1 − z_p = y(1−x_A)`, so both peaks sit at `z → 1`
  and `D(z) → (2α/π)·ln(Q²/m_e²)/[y(1−x_A)]` diverges like `1/y`, while the
  t-peak gets no matching growth (`Y₊ → 2`) and is simultaneously crushed
  because `Q² ≥ 20` forces `x·y ≥ 5.0e−3` — low `y` means **high `x`** — and
  the t-peak's own elastic vertex starts at `t_min = M_A²x_A²/(1−x_A) ∝ x²`
  (0.63 GeV² at `x = 0.79`, where the nucleon dipole `G_D² = (1+t/0.71)^−4` is
  already down to 0.079 and still falling like `1/t⁴`). Along
  `Q² = 23.88 GeV²` the ratio is **1.0014 / 1.061 / 3.98 / 2244** at
  `x = 0.01 / 0.10 / 0.30 / 0.72`, tracking `D(z_s) = 0.082 / 1.36 / 4.39 /
  11.5`: it leaves 1 % exactly where `D(z_s)` passes 1, i.e. where one
  emission stops being the right expansion. **That corner is a breakdown of
  `t-peak+ll`, not evidence that `t-peak` is low by ×6444**, and the whole
  excess there is *quasi*-elastic — ⁶Li's coherent form factor is dead at
  `Q′² ≈ Q² ≈ 28 GeV²` and `ll_peaks_spin1` returns exactly zero.

  **Both tails are tiny there** — at the worst cell `w_tail − 1` is
  4.92e−08 (`t-peak`) and 3.17e−04 (`t-peak+ll`) — which is why the
  event-weighted mean survives. What does **not** survive is the **tensor
  fraction of the tail**, which collapses **−1.41817e−08 → −2.20067e−12
  (×1.55e−04)** at that cell, below the smallest entry of §B2.1's tensor
  table (2.59e−04, and that one is at `Q² = 3`) and **inside** the
  `Q² ≥ 20 GeV²` window. The *absolute* tensor term does not move (both
  models give −6.975e−16); the denominator does.

  It is false outright elsewhere: at
  the HERMES deuteron point (`x = 0.012`, `y = 0.85`, `Q² = 0.53`) the t-peak
  is only **23 %** of the leading-log total — low by a factor **4.36** — and at
  the low-`Q²` corner of the generator window (`x = 0.01`, `y = 0.1`,
  `Q² ≈ 4 GeV²`) the quasi-elastic s+p is **3.35×** the t-peak unsuppressed
  and **7.09×** at the shipped Pauli `k_F` (the s-peak **alone** is 2.33× /
  4.94×). Do not port the `t-peak` default to fixed-target kinematics without
  the s-/p-peaks.
* At θ_S ≠ 0 the tail is the **φ-integrated** one applied to a density that
  carries cos φ′ and cos 2φ′ modulations: it dilutes the φ-*averaged* rate
  correctly and the φ-*differential* rate only on average. Bin in φ at a
  non-longitudinal axis and `rc_tail` is a bin-integrated correction, not a
  per-φ one.

The run banner prints the band at two x values and, **only on a channel where
the tail actually applies**, the form factor in use with its provenance, the
clipped tail-**node** fraction globally and per y-band, and the t-peak caveat
above; on a channel where the tail is off it prints the exclusion reason once
instead of three times. After generation it prints the clipped **event**
counts for the tail and the band — a different quantity from the node
fractions, since nodes are not event-weighted — and both land in the npz as
`meta["rc_clipped_tail_event_fraction"]` / `..._band_event_fraction"]`, with
the per-event bits on `Event.rc_clipped` (`1` = tail, `2` = band). At the
current defaults **nothing clips on the tail** (the pre-2026-09-03 numbers,
1.71 % of nodes and 21.95 % of the `y > 0.9` ones, were an artefact of the
factor-6 per-nucleon error) and 1.2–2.6 % of tagged events clip on the band.
On `--rc-tail-model t-peak+ll` the **nodes** do clip again — **0.36 %**
globally and **4.62 %** of the `y > 0.9` band on the 101 × 77 CLI grid — which
is the same `y → 1` edge as everything else in this bullet; the **event**
count stays **0 / 200 000**, because no accepted cell centre lands in those
nodes. On `--rc-tail-model polrad-full` the nodes clip **0.386 %** globally,
and in a **different** corner — **0.456 %** of `y < 0.5` and **0.124 %** of
`0.5 ≤ y ≤ 0.9`, and **nothing** at `y > 0.9`. Six accepted cells reach the
`tail_max` ceiling there, all at `x = 0.955` with `y` = 0.049–0.193 (**4.96e−07
of the cross section**, and the event count is again **0 / 200 000**): at
`x → 1` the DIS Born is negligible while the elastic tail is not, so `w_tail`
is genuinely ~4500 there and the t-peak, whose own vertex sits above
`t_min ∝ x²`, returns 6.6e−21 — it is not a lower bound with an error there,
it is zero where the answer is everything.
`Pipeline::rc_model()` (Python: `p.rc_model`) exposes `delta(x)`,
`tail_ratio_at(x, q2, q_n)`, `tail_sigma_at(x, q2)` (**five** wide since the
2026-09-03 run: `u, t, qe, u_sp, qe_sp`, the last two zero unless
`t-peak+ll` — and zero under `polrad-full` too, for a **different** reason:
Eq. (18) puts all three peaks in ONE integral, so `u`/`qe` already contain the
s+p and there is no decomposition to read), `ff_provenance`,
`clipped_fraction_by_y` and the **five** raw tail tables for plotting.

Cost: **0.106 s** at setup (the η_A quadratures over a 101 × 77 node grid;
**1.29 s** on `--rc-c0-shape vmc-ft`, which evaluates a 99-term j₀ sum at every
η node below its q = 3 fm⁻¹ cut instead of a closed form; **9.84 s** on
`--rc-tail-model polrad-full`, which replaces each of them by a four-panel
tanh-sinh τ_A quadrature in `long double`) and
**≈ 14 %** of the inclusive event rate (575 k → 495 k ev/s single core on this
machine — the same measurement `OPEN_ITEMS_SOLUTIONS.md` §9 and
`phase_C_numbers.md` §8.3 quote; the *ratio* is the number to carry, the
absolute rates are machine-dependent).

## 7c. Knob provenance — what this run READ, and what it did not

Every run — every channel, every plan, the all-default one included — now
answers one question in one place: **which of its knobs could have affected
the file it just wrote?** `Pipeline.knob_provenance(context)` returns one row
per user-settable knob (**69** on the shipped default `lipolgen-run`, and 69
from `lg.run()` since 2026-09-16, when that entry point started passing
`pzz_mode` — it reported 68 before, missing exactly that row; 65 with an
EMPTY `KnobRunContext`, which is what a bare `Pipeline.knob_provenance()` in
a notebook gets, because the four run-plan rows and `pzz_mode` are things the
core cannot see):

| field | what it is |
|---|---|
| `name` | the `meta` key stem / library field (`pol_sf`, `rc_fq_scale`, …) |
| `flag` | the CLI switch (`--pol-sf`), empty for an API-only knob |
| `value` | this run's own value, as recorded |
| `status` | `read`, `not-read` or `refused` — see below |
| `reason` | one sentence about **this** run, never empty |
| `label` | the short scope clause of a `not-read` reason — the string a `meta` key carries in place of the value |
| `at_default` | whether `value` is the shipped default |

**The three statuses.**

* **`read`** — the run consults it: some quantity it computes is a function of
  this knob, so another value would in general give another file.
* **`not-read`** — nothing the run computes is a function of it, *and a
  different value is accepted*. The `meta` key then carries the `label`
  (`not read on channel coherent-6Li`, `not read by plan tensor-thirds`,
  `not read at rc = off`, …) instead of a value that would mislead.
  **The three route knobs are the one place where the label carries the value
  too** (`--optics`, `n_sigma`, `pot_config`, since 2026-09-05): on a channel
  that writes a far-forward fragment the classifier IS consulted — every
  event's `route` is priced at THIS envelope — and what is absent is only the
  SENSITIVITY of the labels to it on this sample, so the label reads
  `not read on this run at 10x100 high-acceptance (the classifier IS consulted
  here; this run's route labels are insensitive to yr-high-divergence)` and the
  envelope name survives into `meta["optics"]` and the banner's header line.
  On the inclusive channel, where `route_of` returns `Route::Lost` before it
  looks at an envelope, the label stays the bare scope clause
  `not read on channel inclusive`.
* **`refused`** — the **axis** is closed on this run: `PipelineConfig.validate()`
  (or the `Pipeline` constructor) throws on any value but the one shown, so a
  bare value cannot mislead and is kept.

**Where it shows up.**

1. **The npz / HFS `meta`** gains one key, `knob_provenance`, a mapping of
   `name → {value, status, reason, flag, label, at_default}`. It is written on
   every run and survives the npz round trip (the `meta` is JSON). Every
   pre-existing key keeps its name, its type and its meaning; five of them —
   `pol_sf`, `rc_scope`, `coherent_t2`, `pom_set`, `pom_rescale` — now carry
   the **label** where the knob did not run, exactly as `b1_model` has carried
   `none (spin 3/2: no rank-2 input)` since 2026-09-04. `pom_set` stays the
   int it always was wherever it ran.
2. **The run banner** prints one `KNOB PROVENANCE` block: the read knobs
   (values for the ones off their defaults), then every **not-read knob that
   is set away from its default with its reason spelled out** — those are the
   ones that would mislead — then the not-read defaults grouped by their
   label, then the refused axes. Reasons for the rest are in the file.

`lipolgen-run --events 200 --seed 5 --pol-sf nnpdfpol`, verbatim except for
the two `…` elisions:

```
  KNOB PROVENANCE -- what this run READ and what it did not (68 knobs;
     meta["knob_provenance"] carries the whole table with every reason)
     read (16), of which set away from the default: seed = 5, events =
       200, pz = 0.7, pzz = 0.6
     NOT READ, and set away from the default -- a value here would
       mislead, so the file records the LABEL and the reason:
       pe = 0.7  (--pe)
         not read by plan tensor-thirds: every category is built at lam_e = 0
         (bookkeeping.cpp: the three tensor plans hard-code an unpolarised
         beam), so P_e multiplies nothing -- the per-event `pe` column records
         0
       pol_sf = nnpdfpol  (--pol-sf)
         not read under this run's fill: no category carries lam_e * P_e != 0,
         and lam_e * P_e is the only thing g1 is multiplied by
         (InclusiveKernel::amplitudes adds helicity * (m/J) * cos(theta_S) *
         A_par and nothing else reads a PolSF).  The three TENSOR plans build
         every category at lam_e = 0, pe = 0 (bookkeeping.cpp: ...) ...
     not read, at their defaults (44):
       not read at fsi = off: fsi_sigma_mb
       not read at rc = off: rc_delta_low_x, rc_delta_high_x, rc_x_low, ...
         ... rc_tail_model, rc_sp_tensor_scale, rc_with_qe_tail, rc_with_tail,
         rc_n_eta, rc_tail_max, rc_m_lepton
       not read by plan tensor-thirds: pzz_mode
       not read in FIXED-COUNT mode (events = 200): lumi_pb
       not read in fixed-count mode: poisson, apply_optics_lumi_fraction
       not read on channel inclusive: optics, n_sigma, pot_config,
         cluster_beta, p_d, triton_sf, tier, inclusive_b1, coherent_t_max, ...
       not read without --hadronize: coherent_t2, pom_set, pom_rescale
     refused axes (7), where validate() throws on any other value:
       b1_band_scale, b1_alpha_d_dwave_weight, b1_unpol, r_source,
       cluster_wave, cluster_vmc_mc_sigma, fsi
```

*(66 knobs and seven refused axes until 2026-09-06, when `--r-source` added
an eighth and `--pzz-mode` a 68th knob — labelled, not refused, under the
three tensor plans; re-measured from the command above. The eighth refused
axis was `rc_m_lepton`, and it came back OUT the same day: `PolradFull` was
implemented, `m_lepton` is READ on it, and at `--rc off` the row is
`not-read` like every other rc knob — validate()'s whole rc block sits inside
`if (rc != PipelineRc::Off)`, so `refused` was never true there. That is the
ONE `meta` change an `--rc off` run sees from
`docs/open_items/run_2026-09-06/phase_B_numbers.md` §B2; the 48 numeric
columns are byte-identical.)*

Read that block against the five rounds: `pol_sf` and `pe` are round 2, on the
**default plan**; `optics` / `n_sigma` / `pot_config` are not read on the
inclusive channel because `route_of` returns `Route::Lost` before it looks at
an envelope; `cluster_beta`, `p_d`, `triton_sf`, `inclusive_b1` and
`coherent_t_max` are round 5, the ones that used to be recorded nowhere; and
`coherent_t2` / `pom_set` / `pom_rescale` are round 3.

3. **`python/tests/test_knob_provenance.py`** rebuilds the (spec × knob)
   matrix — 12 (isotope, channel, plan) specs, 81 knob variants, **627 cells**
   (re-measured 2026-09-15; 76 / 592 was the count after Phase A)
   (re-measured 2026-09-06, after `--r-source` added three and `--pzz-mode`
   two):
   every knob on the 7 (isotope, channel) combinations under one plan each,
   plus the plan axis on inclusive-⁶Li and tagged-⁷Li-alpha for the 11
   plan-sensitive knobs, which is *not* the full (channel × plan) product —
   and asserts the table against the **output hash** of two small runs: moved → `read`, did not move → `not-read`, refused → `refused`. It
   also fails if a table row has no matrix entry and no excuse, so a knob
   cannot be added without classifying it.

**This block replaced the per-knob reach sentences**, which were true of the
knobs somebody had thought to write one for and silent about the rest. Two are
gone: the SF block's `--unpol-sf reaches EVERY kernel this run builds …;
--pol-sf reaches the inclusive and the tagged kernels`, which was **false on
the CLI's own default plan** (see §2b), and the `--pol-sf DID NOT RUN HERE`
block, which covered one channel of one knob. The PomSet **band** stays where
it was — that is a physics number, not a reach claim.

**Refuse or label — the criterion, and it is stated once** in
`KnobProvenance`'s header block (`include/lipolgen/pipeline.hpp`). A knob this
run does not read is REFUSED when its value would name a *variation of a piece
that did not run* — a scale, a band edge or a shape on a term the run computes
as identically 1 — because such a value claims a systematic was PRICED, and no
label makes a priced systematic un-priced. It is LABELLED when it names a
*backend, an axis or a member of a family that a channel-, plan- or set-scan
sets uniformly*, because refusing one cell of such a scan costs more than it
buys. A rule that depends on the run PLAN can only be labelled, because
`validate()` has no plan.

**What that criterion changed here** (nothing at any default; every generated
array is bit for bit):

* `--rc-fq-scale`, `--rc-tail-tensor-scale`, `--rc-c0-shape`,
  `--rc-qe-suppression`, `--rc-qe-tensor-scale`, `--rc-sp-tensor-scale`,
  `--rc-tail-model` and the four
  API-only tail knobs (`with_qe_tail`, `n_eta`, `tail_max`, `m_lepton` — the
  last of which is **read** on `--rc-tail-model polrad-full` since 2026-09-06
  and refused on the two t-peak models) are **refused** on the three tagged
  channels, where
  `rc_tail == 1` by construction — each of them was measured bit-identical to
  the `--rc tensor-band` baseline there — and every rc sub-knob is refused on
  the coherent channel, where `RcModel::applies()` is false. `--rc
  tensor-band` itself stays accepted on every channel and still prints why it
  prices nothing.
* `rc_options.scope` is **labelled**, not refused, when the fill sits at
  θ_S = 0: `tensor-all` differs from `tensor-rate` only in the cos 2φ
  amplitude, which carries a sin²(θ_S). It is bit-identical to `tensor-rate`
  under `tensor-thirds` and does move under `transverse-tensor` /
  `tensor-flip`.
* `--cluster-wave`, `--triton-sf`, `--inclusive-b1`, `--cluster-beta`,
  `--p-d`, `--fsi-sigma-mb`, `--coherent-t-max` and `--x-max` were accepted
  and recorded **nowhere**; they now have a row on every channel, so a VMC and
  a Hulthén tagged file are no longer indistinguishable in `meta`.
* `--pom-set` / `--pom-rescale` / `--coherent-t2` are labelled off the
  coherent channel and without `--hadronize`.

The full measured reach table, channel by channel and plan by plan, is
`docs/open_items/run_2026-09-03/phase_D_numbers.md` §D6.3–§D6.5.

## 8. Command-line generators

```bash
./build/generate_inclusive --isotope 6Li --config 1 --plan azz --events 200000
./build/generate_tagged    --channel 6Li-alpha --optics tagging --events 400000
./build/generate_tagged    --channel 7Li-alpha --plan apar --events 300000
./build/generate_coherent  --config 1 --slope 50 --f0 0.04 --events 300000
```

Each prints σ per spin category, the tag fraction at every optics of the menu,
the conservation residual and the throughput, and writes HepMC3 (`--out ""` to
skip). Measured single-core throughput on this machine (200 k events, 6Li mid
configuration), with `Pipeline::event` reconstructing in place since
2026-08-30:

| channel | streaming (`for_each`, `for_each_range`) | `generate_range`, chunk 4096 | chunk 1024 |
|---|---|---|---|
| inclusive | **2.83 M ev/s** (was 2.67) | 2.20 (was 2.07) | 2.62 (was 2.36) |
| tagged ⁶Li α at T1 | **1.02 M ev/s** (was 0.93) | 1.01 (was 0.92) | — |
| coherent | **2.32 M ev/s** (was 2.21) | 2.33 (was 2.03) | — |

A `generate_range` buffer costs cache, not allocations, so prefer
`for_each_range(sink, first, last)` when the records are consumed and
dropped. HepMC3 output is the bottleneck when it is on (~18 k ev/s), and the
T2 (PYTHIA) tier runs at ~40 k ev/s. Threading helps only when the per-event
work is heavy; on bare T0/T1 the events are already ~1 µs.

## 9. Polarized ⁶Li configurations for coherent-diffraction codes

`include/lipolgen/cluster_config.hpp` is **opt-in and inert**: nothing in the
generator calls it, it adds no default and changes no existing number. It
produces *nucleon-position configurations* of a polarized ⁶Li — N × 6 × (x, y,
z in fm, isospin) — the input a Good–Walker dipole-model code
(`hejajama/subnucleondiffraction`, arXiv:2408.13213) averages the coherent
amplitude over. It is **not** a diffractive amplitude and not a cross
section; `CoherentScenario` keeps its scenario numbers unchanged, and this
module is the tensor cos 2Φ half of open item 11 (`OPEN_ITEMS_SOLUTIONS.md`
§11 has the unpolarized eSTARlight half and the full derivation).

**Physics, in one paragraph.** ⁶Li(1⁺) is a rigid α(0⁺) core plus a deuteron
with relative L = 0, 2 coupled to S = 1 — the same recoupling `tagged.hpp`
already implements, read in r space instead of k space. For each ion
substate m the deuteron projection m_S is drawn first, then (R, cos θ_R)
from the **m_S-conditioned** table |A_{m_S}(m; R, c)|² (`build_amp2`'s own
discipline, not the m_S-summed density — `rho_alpha_d_summed` exists for
plots only and must not be sampled from), then the p–n pair from AV18 u/w
for that same m_S; the α's four nucleons come from the ANL VMC ⁴He one-body
density, drawn independently and then recentred so Σ s⃗_i = 0 exactly (worst
component 2.2e-15 fm over 20 000 configurations). That is the same diagonal
truncation `tagged.hpp` makes and documents: only the *relative* azimuth of
R̂ and r̂ is lost, and neither ⟨r²⟩ nor the quadrupole Q is affected by it.
The three m-state densities differ only through the Clebsch–Gordan
recoupling C_L(m, m_S) — the same radial tables feed all of m = +1, 0, −1.
Throughput is 0.54 µs/config in-process and 1.31 µs/config through the CLI
(1.06 s wall for 10⁵ configurations including a 37 MB write and its md5),
both against a < 5 µs/config target.

**CLI and Python entry.**

```bash
./build/lipolgen-configs --m +1 --n 100000 --theta-s 1.5707963 --phi-s 0 \
    --seed 20260902 --out li6_m+1.dat --moments-json li6_m+1.moments.json
./build/lipolgen-configs --m unpolarized --n 100000 --out li6_unpol.dat
```

`--m {+1,0,-1,unpolarized}` (`1` also accepted for `+1`),
`--alpha-source {vmc,gaussian}`, `--alpha-d-source {fit-rescaled,fit-raw,overlap-raw}`,
`--alpha-d-scale`, `--match-li6-radius`, `--quadrupole-target`,
`--min-nn-separation`, `--format {he3,annotated}`. Installed as the console
script `lipolgen-configs`; `python -m lipolgen.configs` is the same thing.
From Python:

```python
import lipolgen as lg
s = lg.ClusterConfigSampler()                 # transverse axis by default
st = s.sample_set(100000, m_ion=1, seed=20260902)
st.positions.shape                            # (100000, 6, 3) float64, fm
st.isospin.shape                              # (100000, 6) int64
s.q_matter_analytic_fm2(1), s.eps_b0_equivalent()
```

**The m-state and axis conventions.** `m_ion` is the ⁶Li substate the
configuration was drawn for, {+1, 0, −1}, or the interleaved unpolarized mix
(`sample_set_unpolarized`, reported as m = −2, equal thirds to within one
count). The quantization axis is `(theta_s, phi_s)` in the ion rest frame,
applied as R_z(φ_s) R_y(θ_s) — the identical convention `spin.hpp` /
`tagged.hpp` use for `boost_spectator`. The default is **transverse**,
θ_s = π/2, because the cos 2Φ signal needs it: with the axis along the beam
the projected density is azimuthally symmetric and a₂ ≡ 0 for every m — not
an error, `delta_perp_analytic_fm2` returns exactly 0 for a longitudinal
axis by construction.

**The output format and where it is consumed.** One line per configuration,
byte-compatible with `he3.dat`:

    x1 y1 z1 ... x6 y6 z6   t1 ... t6   m

18 coordinates in fm (ion rest frame, c.m. at the origin, the quantization
axis already applied), six isospins (+1 p, −1 n; nucleons 1–4 are the α core,
5 the proton and 6 the neutron), then the substate. The consumer —
`Nucleons::InitializeTarget` in `subnucleondiffraction` — reads only the
first 3A fields; the rest are invisible to it, exactly as `he3.dat`'s own
trailing fields are. **No comments** — the upstream reader does a bare
`ss >> x` per field and a `#` line would silently produce a wrong
configuration. The metadata therefore lives in the mandatory sidecar
`<out>.meta.json`: the axis, the substate, seed/run and RNG stream, all 13
option fields, every input table with its size, md5 and printed
normalization, the analytic and sampled moments, `quadrupole_dial_s`, and
the quadrupole band. Each input entry also carries `"used"`, so a run with
`--alpha-source gaussian` or `--alpha-d-source overlap-raw` says which of
the four tables it actually read; the caveat block's overshoot factor is
computed from that run's own band, not restated. With `--format annotated`
the same JSON is *also* copied into the `.dat` as a `#` header block — but
that copy is written by C++ before `configs.py` runs, so its `md5`/`git`
stay `null`: **the sidecar is the authoritative record**. A worked example (100 configurations, m = +1, seed 1)
is `docs/open_items/run_2026-09-02/example_li6_m1_configs.dat` plus its
`.meta.json`. Reading the upstream code as it ships today requires the
`-configfile/-configid` generalization of its `A == 3` branch — see the
collaboration ask in `OPEN_ITEMS_SOLUTIONS.md` §11.

**The quadrupole dial, and its floor.** `--quadrupole-target` /
`ClusterConfigOptions::quadrupole_target_fm2` rescales the α–d D-wave
amplitude by a root s of a quadratic (design (G9)) so that the geometry's
Q_charge lands on the requested value — **linear in the D amplitude through
the S–D interference term (95 % of Q), not `sqrt(target/model)`**, which
misses by 7×. This is a **deformation dial, not a wave function**: it exists
so a downstream consumer can ask for the measured tensor moment without
believing the α+d model's own value. **η is EXACTLY linear in this dial**
(η(s) = s·η(1), to 1e−13), because the dial scales R₂ in place and
`asymptotic_ds_ratio()` divides the two waves afterwards — so reading η back
off a dialled configuration recovers the dial, not the wave function. The reachable range for the default
source is s ∈ [0, 1] mapping to Q_charge ∈ **[−0.615, +0.270] fm²** — at
s = 0 the α–d interference term vanishes and only the deuteron's own
+0.270 fm² survives, at s = 1 the wave functions are unmodified. `validate()`
**throws** outside that band. The root that reproduces the measured
Q(⁶Li) = −0.0818 fm² (`LI6_QUADRUPOLE_FM2`) is **s = 0.4032**, closing to
1e-9, and it drags P_D(α–d) down to 3.3 × 10⁻³ from its natural 0.02011 — a
reminder that "match the quadrupole" and "keep the natural D-state
probability" are not simultaneously satisfiable in this model. The writer
stamps `quadrupole_dial_s` beside the band whether or not the dial is used.

**The numbers** (default `FitRescaled` source; full table and every source
variant in `docs/open_items/run_2026-09-02/phase_G_numbers.md`):

| quantity | value |
|---|---|
| P_D(α–d) | 0.02011 |
| ⟨R²⟩_αd | 16.9633 fm² (rms 4.1187 fm) |
| 𝒬[R₀,R₂] | −1.3204 fm² (interference −1.2573, pure-D −0.0631) |
| ⟨r²⟩(⁶Li) | 6.4447 fm² → r_rms 2.5386 fm |
| Q_matter(±1) / Q_matter(0) | −1.2309 fm² / +2.4618 fm² |
| Q_charge(+1) (this geometry) | −0.6154 fm² |
| δ⊥ per nucleon (transverse axis) | −0.1026 fm² |
| eps_b0 equivalent (B = 52.04 GeV⁻²) | −0.0506 (−0.0527 at `CoherentScenario::slope_b` = 50 — see the ΔB note below) |
| a₂(±1) at \|t\| = 0.3 GeV² | +0.1976 |
| asymptotic η (D/S, Whittaker-divided) | −0.0482 (measured `LI6_ETA_DS_GK` = −0.025 ± 0.006 ± 0.010; `OverlapRaw` gives −0.0538 ± 0.0021 stat) |

**The caveats, which the CLI prints on every run and the sidecar stamps.**
The α+d truncation reproduces the ⁶Li point radius to **3 %** — r_rms
2.5386 fm against the measured 2.4655 fm (`LI6_R2_POINT_FM2`), and 4 % above
the VMC `li6.density` value 2.4433 fm (`LI6_R_POINT_VMC_FM`), which is the
number `match_li6_radius` targets — but **overshoots Q(⁶Li) by a factor
≈ 7.5**:
Q_charge = −0.615 fm² (model range −0.615…−0.730 across the three α–d
sources) against the measured −0.0818 fm² and GFMC AV18+IL7's −0.20(6) fm².
Always quote `quadrupole_band_fm2()`, never one number, and never derive a
published tensor input from these wave functions —
`docs/OPEN_ITEMS_SOLUTIONS.md` §11's rule stands. The asymptotic D/S ratio η
says the excess is a real but *moderate* ≈ 2× effect, not the 5–15× a naive
R₂/R₀ ratio would suggest (§2.1 of `phase_G_numbers.md`).

That 7.5× is **two factors, not three**: **3.3165** from the model to the dial
setting that matches the *measured* η, × **2.2686** from there to the
measurement. The missing ≈15 % non-α+d component (1/S_αd = 1.1706) is **not** a
third factor — both waves are divided by √S_αd before the moments are taken, so
it is already inside the −0.615. **Quote the band, not the leg**: the measured
η carries ±0.011662, and mapped through the (exactly linear) dial that spans
model Q from −0.4005 fm² to **+0.0298 fm² — through zero**. So the 3.32× is
*3.32× (1 σ: 1.54× … sign change)*, and η cannot separate the measured
−0.0818 fm² (+0.48 σ) from GFMC's −0.20 (−0.07 σ).
`docs/open_items/run_2026-09-03/phase_C_numbers.md` §C1 has the derivation;
`LI6_ETA_DS_GK` / `_STAT` / `_SYST` are the single home of the measurement.

**Driving the one dial from η: `quadrupole_for_eta`** (open item C5.1,
2026-09-04). η has deliberately **no dial of its own** — it is exactly linear
in `quadrupole_dial_s()`, so an `eta_target` option would be a second name for
the same one-parameter family, which `docs/CONVENTIONS.md` forbids and which
could be set to contradict `quadrupole_target_fm2` in one options object. What
was missing was the *conversion*, which lived only as a typed number in a
document; it is now a method:

```python
s = lipolgen.ClusterConfigSampler()
q = s.quadrupole_for_eta(lipolgen.LI6_ETA_DS_GK)   # -> -0.1842160147 fm^2
o = lipolgen._lipolgen.ClusterConfigOptions(); o.quadrupole_target_fm2 = q
lipolgen.ClusterConfigSampler(o).asymptotic_ds_ratio()   # -0.025, to 1e-12
```

It uses the sampler's own **pre-dial** moments and its current dial, so it is
exact on an already-dialled sampler, and it **throws** rather than clipping
when the requested η is out of reach. Feed it the GK band and it re-derives
the sign change from the tables: −0.0366619 → **−0.4004752 fm²**, −0.025 →
−0.1842160 fm², −0.0133381 → **+0.0297575 fm², the wrong sign**.

The default α–d source is `FitRescaled` by **author decision**, not because the
smoothed `li6.adr.fit` R₂ node near 1.07 fm is preferred: that node is **real**
in the raw block (3.3 σ) and irrelevant either way (r < 1.5 fm is −0.09 % of
q_int). `OverlapRaw` moves Q_charge by +8.7 % and η by +11.5 %, inside a band
that is already a factor 7.5 wide; `phase_C_numbers.md` §C3 has the cost of
each choice. The α core is an
**uncorrelated** product of one-body densities (`min_nn_separation_fm = 0`
by default). Switching the hard core on (`--min-nn-separation 0.9`, what
arXiv:2605.00454 imposes) **moves ⟨r²⟩ by ≈ +2 %** — 6.573 fm² against the
free 6.445 — because rejection correlates the four s⃗_i and the closed form
⟨s²⟩ = (3/4)⟨v²⟩ holds only for *independent* draws. The constructor
therefore measures ⟨s²⟩ once (10⁵ cores on a fixed stream), stamps
`alpha <s^2> NUMERICAL (hard core)` into the provenance and
`moments.r2_mean_is_approximate` into the sidecar, and the CLI prints a
note; `rho_alpha_recentred()` is *not* corrected and remains the
independent-draw closed form. The exact 5-D joint density with the m_S
angular coherences restored is `ClusterConfigOptions::exact_coherence`,
which **throws** — reserved, not implemented, so the approximation is
visible in the API rather than buried in a comment.
`docs/open_items/run_2026-09-02/design_G_cluster_config.md` has the physics
and `phase_G_numbers.md` the reproduced numbers.
