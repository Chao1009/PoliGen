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
Event ev;  while (p.next(ev)) { ... }             // pull API (cursor; rewind())
std::vector<Event> batch;                          // explicit index range
p.generate_range(1000, 2000, batch);
```

Bookkeeping, all share-invariant in pb:

```cpp
p.sigma_per_category_pb();   // [pb], no lumi_fraction, no Scenario::run_share
p.sigma_pb();                // sum_k lumi_fraction_k * sigma_k
p.lumi_per_category_pb();    // [pb^-1] (luminosity mode)
p.counts();  p.size();       // events per category / in total
```

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
                                          // rank-2 transfer + 1e-2 toy Delta
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
  if (rp_tagged(ev, p.optics(), p.pot_config()))         // Roman-Pot mask
    printf("k=%.3f cos=%+.3f m_S=%+g alpha_s=%.3f pT_s=%.3f\n",
           ev.kin.k, ev.kin.cos_theta_k, ev.spin.m_struck,
           ev.kin.alpha_s, ev.kin.pt_s);
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

Balance: **whole nucleus**, `k + P_ion = k' + p_spec + X`. The spectator is on
shell at its AME2020 mass, the struck cluster takes the remainder
`P_X = P_ion − p_spec` and is off shell, and `X = k + P_X − k'`. Charge closes.

## 4. Coherent ⁶Li

```cpp
PipelineConfig cfg;
cfg.channel = PipelineChannel::CoherentLi6;
cfg.isotope = "6Li";
cfg.n_events = 300000;
cfg.coherent.f0      = 0.04;     // band 0.02 - 0.08     (SCENARIO)
cfg.coherent.slope_b = 50.0;     // band 40 - 60 GeV^-2  (SCENARIO)
cfg.coherent.amp     = 0.01;     // flat cos2phi at P_zz = 1, band 3e-3 - 1e-2
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

It is called once per finished T0 event with the event's own counter-based
stream, after every T0 particle is in place. For a tagged event the hard process
belongs to the **struck cluster**: `Role::StruckCluster` carries the off-shell
`P_X`, `kin.alpha_s` / `kin.pt_s` are the light-front variables in the ion rest
frame, and

```cpp
StruckCluster sc;
struck_cluster_of(ev, *p.tagged_channel(), sc);
sc.p_per_nucleon_eff;   // P_X,z / A_partner -- inject the hard process here
sc.virtuality;          // M_X^2 - m_free^2 < 0
```

rebuilds the whole record. A hadronizer used with `for_each(sink, nthreads)`
must be re-entrant.

## 7. Checks you get for free

```cpp
Vec4   r = momentum_residual(ev);   // zero to rounding, per the channel's balance
Vec4   s = momentum_scale(ev);      // what entered that balance
double q = charge_residual(ev);     // exactly zero
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
skip). Measured single-core T0 throughput on this machine: **3.5 M ev/s**
inclusive, **1.5 M ev/s** tagged, **2.7 M ev/s** coherent; HepMC3 output is the
bottleneck when it is on (~18 k ev/s). Threading helps only when the per-event
work is heavy (the T2 tier); on bare T0 the events are already ~700 ns.
