# LiPolGen — the C++ API in two pages

Everything below is `namespace lipolgen`, C++17, header `<lipolgen/pipeline.hpp>`,
link `LiPolGenCore` (+ `LiPolGenHepMC` if you write HepMC3).

```bash
source env.sh
cmake -S . -B build && cmake --build build -j8 && ./build/lipolgen_tests
```

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
