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
built on, at that same path. Making it relocatable means running
[`auditwheel repair`](https://github.com/pypa/auditwheel) after `pip wheel`,
which vendors HepMC3/PYTHIA8/LHAPDF's shared libraries into the wheel itself
and rewrites the RPATHs to be `$ORIGIN`-relative — standard practice for
compiled-extension PyPI wheels, but note it turns the wheel into a genuine
combined/linked work with all three (GPL-2-or-later, GPL-3.0, GPL-3.0), so
whoever redistributes an `auditwheel`-repaired wheel is redistributing under
GPL-3.0-or-later terms for the combination, regardless of LiPolGen's own
license header (see `docs/OPEN_ITEMS_SOLUTIONS.md` §12–13 and `docs/open_items/engineering.md` §C).

## 1. The shape of a run

A run is **one `PipelineConfig` + one `RunPlan`**.

* `PipelineConfig` says *what physics and how much of it*: isotope, beam
  configuration index (0/1/2 = low/mid/top), channel, acceptance window
  (`Scenario`), far-forward optics, seed, and either `lumi_pb` or `n_events`.
* `RunPlan` (`bookkeeping.hpp`) says *how the luminosity divides between spin
  fills*: `tensor_thirds_plan`, `helicity_flip_plan`, `transverse_tensor_plan`,
  `tensor_flip_plan`, or your own `std::vector<SpinCategory>`.

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
| `cdks` | `B1Model::Cdks` | the same ⁶Li rank-2 transfer on the **other camp** for b₁ᵈ — the digitized Cosyn–Dong–Kumano–Sargsian PRD **95** (2017) 074036 Fig. 4 column. \|b₁\| is **two orders of magnitude smaller below x ≈ 0.1 and comparable above it**: peak \|x·b₁\| 1.67e−4 against Miller's 4.27e−4 over x ∈ [0.02, 0.95] at Q² = 2.5, and at x = 0.3 it is **4× larger** with the opposite sign (∫b₁ dx 6.9e−5 against 9.1e−4). Different sign structure throughout. |
| `li6-convolution` | `B1Model::Li6Convolution` | `Li6ConvolutionB1` (`b1_nuclear.hpp`): the **four-term α–d convolution** of `docs/open_items/run_2026-09-02/design_D_b1_li6.md`. |

Miller (HERMES-like) and CDKS (convolution) are *different camps* for the
deuteron's own b₁ and the library does not adjudicate between them. **Say which
one a plot used.**

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
  species. `7Li` is spin 3/2 and has no rank-2 input here.
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

### Which R each object uses — they are NOT the same

`Li6ConvolutionOptions::r_func` null means **`r_sigma_lt`**, the kernel's own
toy R, so that the ⁶Li backend's F₁ is consistent with `InclusiveKernel`'s.
`DeuteronConvolutionB1::Options::r_func` null means **`r1998`**, the SLAC world
fit, because that is what CDKS used and the A = 2 gate exists to reproduce
*their* figure. So **the R the gate was validated with is not the R the shipped
⁶Li backend runs with.** Measured on the gate (κ = 1, like for like): switching
r1998 → r_sigma_lt moves the second zero from 0.392 to 0.365, i.e. *away* from
the digitized 0.457, and the peak up by 1.6 %. Far inside the 100 % band, but
not cosmetic; whether the ⁶Li default should become `r1998` is an **open
follow-up, not decided in the 2026-09-03 close-out** (`phase_D_gate.md`
checklist item 2; `OPEN_ITEMS_SOLUTIONS.md` §10 item 4).

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

### ⚠ The A = 2 validation gate is NOT fully passed

`li6-convolution` is in the tree, opt-in and tested, but the gate of design §5
**fails its magnitude clause**. Fed the AV18 deuteron u(k), w(k) instead of the
α–d waves, the same kernel must reproduce the digitized CDKS Fig. 4:

| clause | result |
|---|---|
| G0a–G0g (analytic identities), G1a–G1f (unpolarized convolution), layer 2 (α–d sign gate) | **PASS** |
| **G3a** — two sign changes and the large-x peak position | **PASS** at the gate's default (κ = √(1+γ²), ToyF2, `r1998`) — but by a thin margin, below |
| **G3b** — peak magnitude within a factor 2 of the digitized 1.0852e−3 | **FAIL**: 4.7748e−4, ratio **0.440**, a factor **2.27 low** (2.9484e−4, ratio 0.272, a factor 3.68, at CDKS Eq. (17)'s κ = 1 form) |

**G3a's margin, stated.** The low-x zero clears the [0.02, 1.0] counting
window's lower edge by **0.0021**, and with a realistic PDF (CT18NLO) it drops
below the scan floor, so the "exactly two sign changes" clause *fails* there.
The low-x zero is **not** a robust discriminator; the second zero (0.392 → 0.377
against the digitized 0.457) and the peak position (0.755 against 0.766) are.

**Since 2026-09-03 the gate is quoted at CDKS Eq. (21)'s δ-function**
(`DeuteronConvolutionB1::Options::finite_q_delta`, now `true`): Eq. (21) is
their exact definition, y = (E − p_z κ)/M_N with κ = \|q⃗\|/ν, and Eq. (17)'s
(E − p_z)/M_N is the "≃" of their Eq. (18). It is worth ×1.62 at the peak.
`Li6ConvolutionOptions::finite_q_delta` keeps its `false` default, where the
same switch is worth 1 %.

The remaining factor 2.27 is *attributed*, not dangling: a real nucleon PDF
(CT18NLO in place of `ToyF2`) is ×1.67, and with it the ratio is **0.719 —
inside the factor-2 window**. That is not this phase's default to change (the
core must not link LHAPDF, and CDKS's own MSTW2008 LO is not installed), so it
is recorded rather than absorbed. The full checklist is
`docs/open_items/run_2026-09-02/phase_D_gate.md`.

**Consequence: no ⁶Li number from this backend may be published while that
stands.** Open item 10 stays open. The CLI prints the warning on every run.

### The top x cell — why `--x-max 0.95`

Both opt-in backends carry the CDKS camp's b₁ᵈ, which is a **Q² = 2.5
digitization with no Q² evolution**. In the topmost cell of the default window
(x = 0.955) F₁ has fallen far enough that b₁/F₁ reaches **3.3** for `cdks` and
**5.6** for `li6-convolution` — past where the phi-averaged density 1 + w_avg
stays positive, and `InclusiveSampler` refuses the run with *"negative
phi-averaged density"*. Use `--x-max 0.95` (or `cfg.scenario.x_max = 0.95`).
`lipolgen-run` now checks this **before** building the pipeline and says which
flag fixes it, rather than letting that message — which names neither `--x-max`
nor the CDKS table — come out of the sampler three frames down.
`miller`'s b₁ is a ratio model (0.145 at the same point) and does not need it.
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
`docs/OPEN_ITEMS_SOLUTIONS.md` §10. **They are recorded for reproducibility and
regression-catching, not as results**, for the gate reason above.

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
scenario placeholder).  The deuteron control channel is always Hulthen.  Data
files are found at `$LIPOLGEN_DATA_DIR`, else the compiled-in
`${CMAKE_SOURCE_DIR}/data`; the tables are zero past 5 fm⁻¹ = 0.9866 GeV, so
no spectator is drawn beyond that.

What it changes (40 k events, 10 × 99.5 GeV/u; full tables and the derivation
in `docs/open_items/vmc_reconciliation.md`):

| | Hulthen β = 0.30 | VMC AV18 |
|---|---|---|
| ⁶Li ⟨k⟩ / P(k>0.2) / P(k>0.45) | 0.1219 / 0.1476 / 0.0157 | 0.1225 / 0.2410 / 0.0019 |
| ⁷Li ⟨k⟩ / P(k>0.3) / P(k>0.45) | 0.2893 / 0.3527 / 0.1555 | 0.1864 / 0.2111 / 0.0114 |
| ⁶Li P_D | 0.0867 | 0.01935 |
| ⁶Li tag fraction, YR high-acceptance | 0.0249 | 0.0348 |
| ⁶Li tag fraction, tagging optics | 0.2530 | 0.2485 |
| ⁷Li tag fraction, YR high-acceptance | 0.9730 | 0.9981 |
| ⁶Li A_zz^tag at k = 0.20 GeV | +0.845 | +0.452 |

**Roman-Pot tag fractions, Hulthén β = 0.30 vs VMC, at the three configurations
(40000 events, seed 20260829):** `validation/vmc_tag_fractions.py --events
40000 --configs 0,1,2` (5×41-class, 10×100-class, 18×275-class nominal design
points; 18×275 rigidity-caps ⁷Li's ion energy below ⁶Li's, so it is quoted
separately per isotope).

| isotope | optics | 5 × 40.8 GeV/u Hulthén | VMC | 10 × 99.5 GeV/u Hulthén | VMC | 18 × 137.5/117.9 GeV/u Hulthén | VMC |
|---|---|---|---|---|---|---|---|
| ⁶Li | YR high-acceptance | 0.0286 | 0.0365 | 0.0249 | 0.0348 | 0.0266 | 0.0349 |
| ⁶Li | tagging optics | 0.3410 | 0.3175 | 0.2530 | 0.2485 | 0.3115 | 0.2905 |
| ⁷Li | YR high-acceptance | 0.9660 | 0.9981 | 0.9730 | 0.9981 | 0.9787 | 0.9981 |
| ⁷Li | tagging optics | 0.9805 | 0.9993 | 0.9927 | 0.9992 | 0.9941 | 0.9992 |

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
```
```bash
python -m lipolgen.cli --rc tensor-band --events 400000
python -m lipolgen.cli --rc tensor-band --rc-delta-low-x 0.19 --events 400000
```
```python
cfg = lipolgen.make_config(events=400000, rc="tensor-band",
                           rc_delta_low_x=0.30, rc_fq_scale=1.0,
                           rc_qe_suppression=1.0)
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
| every `Tagged*` | full (τ from the cluster density), **clamped** | **≡ 1** | the elastic recoil sits at x_L = 1, inside the 10σ beam envelope, so the tag itself vetoes it |
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
* **On the tagged channels the band is CLAMPED, and it has to be.**
  `τ_tag = 1 − n̄(k,c)/n_M(k,c)` is unbounded: wherever the event's own
  `n_M` is near a node of the M-dependent spectator density — the ⁶Li `M = 0`
  density has them — the ratio blows up (measured **|τ| up to 30.7**).
  Unclamped, a 20 k-event tagged-alpha run published `rc_tensor_hi` down to
  **−1.79** (⁷Li: **−8.35**), i.e. **negative weights in the npz**.
  `RcOptions::band_tau_max` (default 1.0) holds every published edge inside
  `[1 − δ, 1 + δ]`; the clipped fraction is **0.6 % of ⁶Li tagged-alpha
  events and 2.6 % of ⁷Li**, is printed by the run and lands in
  `meta["rc_clipped_band_event_fraction"]` and per event in
  `Event.rc_clipped`. `n_M → 0` is exactly where "rescale the tensor part by
  δ" stops meaning anything, so the clamp is a **choice**, not a fix.

**The mandatory band runs — the 20–40 mb rule of `--fsi-sigma-mb`, applied
here.** Never quote one row alone. Run each knob at both/all its ends and
quote the lo/hi envelope on `A_zz`:

| knob | rows to run | what it prices |
|---|---|---|
| `--rc-delta-low-x` | **0.19 and 0.30** | 0.30 is the conservative end of Gakh–Shekhovtsova's 10–30 %; 0.19 is the residual HERMES actually achieved at its lowest-x bin |
| `--rc-fq-scale` | **0, 1, 2** | ±100 % on the ⁶Li quadrupole form factor. σ^el_T is **quadratic** in it, so this band must be **RUN, never rescaled** from one row — the two edges are not symmetric about the nominal, and at `x = 0.01, Q² = 5` they even bracket a **sign change** of ΔA_zz |
| `--rc-tail-tensor-scale` | **0.5, 1, 2** | the η·F_m² tensor sector, which `--rc-fq-scale` does **not** span |
| `--rc-qe-suppression` | **0, 0.5, 1** | a flat multiplier on the quasi-elastic tail, **on top of** the Pauli suppression below. `rc_tail` is quasi-elastic-**dominated** at every `x ≳ 0.03` (73 % of the tail at `x = 0.1`, **99.9 %** at `x = 0.30`), so 0 is never a small variation |
| `RcOptions::qe_kf_gev` (API) | **0.169 (default), 0.221, 0** | POLRAD Eq. (44)'s `S_E`/`S_M`, the de Forest–Walecka Fermi-gas factor `S(q) = (3/4)(q/k_F) − (q/k_F)³/16` below `q = 2k_F`, exactly as POLRAD's `ffquas` codes it. **On by default** at ⁶Li's measured `k_F` (Moniz *et al.*, PRL **26** (1971) 445). It cuts the QRT to **0.47** at `x = 0.01` and **0.87** at `x = 0.1`; `0` is the unsuppressed edge v0 shipped |

At `x = 0.01, Q² = 5 GeV²` (config 1, `--plan tensor-thirds --pzz 0.6`) the
band is by far the largest of these: half-width `4.4e−04` on `A_zz` at
δ_low = 0.30 (`2.8e−04` at 0.19), against `8.0e−07` from the whole `fq_scale`
band and `3.5e−07` from the entire tail. **The band — the unapplied
lepton-vertex correction — is what `--rc tensor-band` is for; the tail is a
bookkeeping item at EIC energies and the headline at a fixed target.** Every
number is in `docs/OPEN_ITEMS_SOLUTIONS.md` §9.

**Honest flags, to repeat wherever any of this is quoted:**

* `δ_low = 0.30` is the **size of a correction this generator does not
  apply**, taken as a 1σ band — *not* a measured residual. Its source
  (Gakh–Shekhovtsova, [hep-ph/0403262](https://arxiv.org/abs/hep-ph/0403262))
  has **zero INSPIRE citations**.
* Everything cited is **deuteron**. Whether the deuteron's *fractional* RC
  transfers to ⁶Li is untested, and is the largest unquantified assumption
  in the band.
* The ⁶Li form-factor shape parameters are **unfitted starting values**, not
  a fit (the elastic data are not in this repository in machine-readable
  form). The *normalisations* are the **measured** moments — `μ = +0.822047
  μ_N`, `Q = −0.0818 fm²` — and never VMC, whose `Q(⁶Li) = −0.23(9) fm²` is
  3× the measured one.
* **No RC calculation exists for a tagged tensor asymmetry**, so the tagged
  band is a defensible but **uncited extrapolation**; and none exists for any
  φ-dependent tensor observable at any axis, which is why the `Δ` cos 2φ
  sector and the coherent channel carry no band at all.
* The **polarised** quasi-elastic tail is not priced at all (Z.-L. Zhou
  *et al.*, PRL **82** (1999) 687).
* **`rc_tail` is ONE PEAK of the elastic and quasi-elastic tails — the
  t-peak — and its absolute normalisation is not validated against any exact
  tail.** It is a **lower bound on the dilution**. POLRAD §2.1.3 B asserts the
  s- and p-peaks are suppressed for a tail; `tests/test_rc.cpp` T8(c) measures
  that with an independent leading-log construction and finds it true **where
  this generator runs** (> 99 % of the total for ⁶Li at `Q² ≥ 20 GeV²`, because
  the ⁶Li charge form factor is dead by `t ≈ 0.25 GeV²` while the s-peak sits
  at `t ≈ zQ²`) and **false elsewhere**: at the HERMES deuteron point
  (`x = 0.012`, `y = 0.85`, `Q² = 0.53`) the t-peak is only **23 %** of the
  leading-log total — low by a factor 4.4 — and at the low-`Q²` corner of the
  generator window (`Q² ≈ 4 GeV²`) the quasi-elastic s-peak is already 3.3×
  the t-peak. Do not port `rc_tail` to fixed-target kinematics without the
  s-/p-peaks.
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
factor-6 per-nucleon error) and 0.6–2.6 % of tagged events clip on the band.
`Pipeline::rc_model()` (Python: `p.rc_model`) exposes `delta(x)`,
`tail_ratio_at(x, q2, q_n)`, `tail_sigma_at(x, q2)`, `ff_provenance`,
`clipped_fraction_by_y` and the three raw tail tables for plotting.

Cost: **0.106 s** at setup (the η_A quadratures over a 101 × 77 node grid) and
**≈ 14 %** of the inclusive event rate (575 k → 495 k ev/s single core on this
machine — the same measurement `OPEN_ITEMS_SOLUTIONS.md` §9 and
`phase_C_numbers.md` §8.3 quote; the *ratio* is the number to carry, the
absolute rates are machine-dependent).

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
believing the α+d model's own value. The reachable range for the default
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
| eps_b0 equivalent (B = 52.04 GeV⁻²) | −0.0506 |
| a₂(±1) at \|t\| = 0.3 GeV² | +0.1976 |
| asymptotic η (D/S, Whittaker-divided) | −0.0482 (measured −0.025 ± 0.006 ± 0.010) |

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
R₂/R₀ ratio would suggest (§2.1 of `phase_G_numbers.md`). The α core is an
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
