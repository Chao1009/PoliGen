# The T2 chain: Pipeline -> PythiaBridge -> HepMC3 -> ePIC

One page: how to run it, what a T2 record contains, what `abconv` / `npsim`
expect of the file, and what still does not fit.

```bash
source env.sh && cmake -S . -B build -DLIPOLGEN_WITH_PYTHON=OFF \
  && cmake --build build -j8 && ./build/lipolgen_tests
./build/generate_full --channel inclusive  --events 200000
./build/generate_full --channel 6Li-alpha  --events 200000 --plan azz
./build/generate_full --channel 7Li-alpha  --events 200000 --plan apar
./build/generate_full --channel coherent   --events 200000 --plan azz
./build/generate_full --channel inclusive  --events 200000 --no-hadronize   # T0 only
./build/generate_full --channel 6Li-alpha  --events 200000 --tier T0        # the old record
```

`--isotope 6Li|7Li` (inclusive only -- the tagged channels imply their own,
coherent is 6Li only), `--config 0|1|2` (energy point), `--seed`, `--out
FILE` (`""` = no HepMC3), `--tier T0|T1` (tagged final-state fidelity,
default T1), `--no-hadronize` (skip PYTHIA, print T0 bookkeeping only). Each
run prints sigma per spin category, the far-forward tag fraction, the PYTHIA
veto/retry/no-surrogate counters, the observed whole-record conservation
residual, and events/s.

**No hook of any kind is needed for a tagged channel any more.** The plain
`cfg.hadronizer = bridge.hadronize` binding conserves the whole record,
because the pipeline resolves the struck cluster before the bridge sees the
event -- section 1a.

## 1. What a T2 record contains

One `lipolgen::Event` per collision. `particles[0]` = beam e-, `[1]` = beam
ion, then T0/T1 (Pipeline) and T2 (PythiaBridge) particles interleaved by
the order each tier appended them.

| role | status | tier | HepMC3 |
|---|---|---|---|
| `BeamElectron`, `BeamIon` | `Beam` (4) | T0 | incoming at the primary vertex |
| `ScatteredElectron` | `Final` (1) | T0 | outgoing, untouched by the bridge |
| `VirtualPhoton` | `Intermediate` (3) | T0 | documentation only (optional) |
| `Spectator` | `Final` (1) | T0 (tagged) | the tagged cluster; outgoing, **never read or written by the bridge** -- bit-identical before/after `hadronize()` |
| `PartnerSpectator` | `Final` (1) | **T1 (tagged)** | the struck cluster's non-struck nucleon(s) / bound remnant; same "never touched" guarantee |
| `IntactRecoil` | `Final` (1) | T0 (coherent) | outgoing; same "never touched" guarantee -- bit-identical through a coherent `hadronize()` |
| `Pomeron` (pdg 990) | `Intermediate` (3) | T0 (coherent) | documentation: the coherent channel's T2 target, `P_IP = P_ion - P_recoil`, spacelike with mass written `-sqrt|t|`; `hadronize()` reads it FIRST and runs the gamma*-Pomeron tier on it (`docs/PYTHIA_BRIDGE.md` §12) |
| `HadronicX` (pdg 92) | `Final` -> **demoted to `Intermediate` (3) by `hadronize()`** | T0 | documentation only once real hadrons exist, so a status-1 sum is never double-counted |
| `StruckNucleon` | `Intermediate` (3) | T0 (inclusive) / **T1 (tagged)**, or appended by `hadronize()` when the target was implicit | documentation: what nucleon 4-vector/pdg actually went into PYTHIA.  Off shell on a tagged event (it is P_X minus the partners) and carrying a +-1 `pol` label |
| `StruckCluster` | `Intermediate` (3) | T0 (tagged) | documentation; at T1 it equals the struck nucleon plus every partner spectator, exactly |
| `Hadron` | `Final` (1) | T2 | PYTHIA's showered/hadronized final state |

`lipolgen::HepMC3Writer` (`docs/HEPMC3_CONVENTION.md`) writes this as one
`GenEvent`: status-4 beams at a primary vertex with both incoming, every
other particle attached by `mother1`/`mother2` (or the primary vertex if
unset), spin/kinematics as named `GenEvent` attributes (`spin_J`, `spin_M`,
`P_z`, `P_zz`, `dis_x`, `dis_Q2`, `spectator_k`, `t`, `channel`, `run`,
`bunch`, ...), and a `GenCrossSection` in pb. `Event::total_final()` /
`total_charge_final()` sum exactly the `Status::Final` rows of the table
above -- that is the identity every conservation check in `test_t2.cpp` is
against.

## 1a. The tagged chain at T1 (the struck cluster is resolved)

`PipelineConfig::tier` is `Tier::T1` by default on the three tagged
channels. Before the T2 hook runs, `ClusterBreakup` (`breakup.hpp`) turns the
off-shell struck cluster `P_X = P_ion - p_spec` into

    P_X  ->  p_N,struck (OFF shell, Role::StruckNucleon, status 3)
           + partner spectator(s) (ON shell, Role::PartnerSpectator, status 1)

| channel | struck cluster | what T1 writes | partners/event |
|---|---|---|---|
| ⁶Li α tag | embedded d | struck p or n (`F2p : F2n`) + the other nucleon | 1 |
| ⁷Li α tag | quasi-free t | struck n + d, **or** struck p + two n | 1 or 2 |
| d control | already a nucleon | the same nucleon, relabelled `StruckNucleon` | 0 |

and `X = k + p_N,struck - k'` becomes the per-nucleon remainder, so

    k + P_ion = k' + p_spec + sum(p_partner) + X            (exact)

with **no caller-side hook**. The bridge's first-priority branch takes
`Role::StruckNucleon` verbatim; its `Role::StruckCluster` branch is
deprecated, warns once and counts into
`PythiaBridgeStats::n_cluster_fallback`.

Measured, 300 events per channel at 10 × 99.5 GeV/u with the plain binding:

| channel | hadronized | worst 4p residual (rel) | worst charge | no-surrogate |
|---|---|---|---|---|
| ⁶Li α | 300/300 | 1.4 × 10⁻¹³ | 0 | 0 |
| ⁷Li α | 300/300 | 5.9 × 10⁻¹⁴ | 0 | 0 |
| d control | 300/300 | 2.2 × 10⁻¹³ | 0 | 0 |

The **"no surrogate" tail is gone** (it was ~5 % for ⁶Li and ~12 % for ⁷Li
with the whole-cluster hook, and the reason is mechanical: the PYTHIA-side
surrogate now has to reach one nucleon's momentum inside its `headroom`
instead of a whole cluster's). `generate_full --tier T0` still reproduces the
old record and prints its honest residual: 0.186 relative, charge wrong by 1.

Two identities move at T1 and one does not:

* `hfs_sigma_empz_exact(k, p_N, k')` -- exact for any target -- still holds
  at 1.5 × 10⁻¹³ GeV.
* the closed-form collinear `hfs_sigma_empz_truth` is now approximate at the
  **1.4 %** level (⁶Li α; 0.7 % for ⁷Li and the d control), because the
  struck nucleon carries the full internal Fermi p_T where the T0 cluster
  carried only the spectator's recoil. That is the formula's assumption
  failing, not the record.
* `Sum p_T,hadrons = p_N,T - p_T,e'` is unchanged, 1.4 × 10⁻¹² GeV.

Mean struck-nucleon virtuality `p_N² - M_N²` over 50 k events:
**-0.057 GeV²** (⁶Li α), **-0.124** (⁷Li α), **-0.035** (d control) -- a
Fermi-motion scale, not a nuclear-mass one. The timelike-X rejection is
redone on the (smaller) T1 X and fires on 0.024 % / 0.022 % / 0.004 % of
draws respectively.

## 2. What `abconv` / `npsim` expect

From `PolarizedLithiumSim/plans/05_doubly_polarized_generator.md` step 5.D
and `tools/fullsim/README.md` (read-only references for this file):

- **HepMC3 Asciiv3**, not HepMC2: DD4hep's `npsim` reads a `.hepmc` file
  with HepMC3's `ReaderAscii`, so a HepMC2 `IO_GenEvent` file is an
  immediate EOF. `HepMC3Writer` only implements `Asciiv3` (`HepMC2Ascii`
  throws at construction) -- already the only option this library offers.
- **Status-4 beams.** `abconv` (the HepMC3 -> Geant4-primary-generator
  bridge) requires the two incoming beam particles at status 4, which is
  what `HepMC3Writer` writes `Role::BeamElectron`/`BeamIon` as regardless
  of `Particle::status` on the T0 record (both are already `Status::Beam`
  by construction).
- **10-digit ion PDG codes** for nuclei (`1000030060` = 6Li, `1000020040`
  = alpha, ...): both the beam ion and any `Spectator`/`IntactRecoil`/
  `PartnerSpectator` fragment already carry these (`event.cpp`,
  `spectator.cpp`, `coherent.cpp`). **`npsim` itself cannot shoot a bare
  ion via `--gun.particle`** ("Bad particle type" -- generic ions are not
  pre-instantiated in the G4 particle table); it has to come from the
  HepMC3 file, which is exactly this chain's output.
- **No bare vertex.** A vertex with nothing incoming (e.g. a single-particle
  gun event written by hand) is dropped by HepMC3's own writer; every
  vertex this library writes has both beams incoming at the primary vertex,
  so this never arises here.
- **Spin labels have no existing HepMC3 convention** (plans/04 #17): this
  library's `spin_J`/`spin_M`/`P_z`/`P_zz`/`spin_axis_theta`/`phi` event
  attributes are a *proposed* convention, not yet an ePIC-MC-group standard
  -- consumers downstream of `npsim`/`eicrecon` that want the spin label
  back need to read the HepMC3 file directly (`abconv`/`npsim` do not carry
  event attributes through to `edm4hep`).
- **Smoke gate** (plans/05 5.D / 5.4): "HepMC3 -> abconv -> npsim,
  event-by-event 4-momentum/charge, 100-event smoke passes" -- this
  library's own event-by-event checks (`test_t2.cpp`) are the upstream half
  of that gate; the `abconv`/`npsim` half is unrun here (no container in
  this environment) and stays an open item below.

## 3. Open items

1. **`abconv` -> `npsim` -> EICrecon smoke (plans/05 5.D) is not run from
   this repository.** `test_t2.cpp` and `generate_full` establish that the
   HepMC3 file is well-formed and event-by-event conserving up to the
   PythiaBridge tier; the 100-event `abconv`/`npsim` pass itself needs the
   `eic_xl`/`jug_xl` containers of `tools/fullsim/README.md`, unavailable
   here.
2. **CLOSED (2026-08-30) -- tagged conservation.** The complete fix named
   in the previous revision of this file ("the cluster's non-struck
   nucleon(s) should be emitted as `Role::PartnerSpectator`") is what the T1
   tier does; section 1a has the measurements. `set_nucleon_in_cluster` and
   the whole-cluster workaround are gone from `generate_full` and
   `test_t2.cpp`, and the bridge's cluster branch is deprecated. What remains
   OPEN inside T1 is physics, not bookkeeping:
   * **FSI.** No final-state interaction of any fragment with any other, no
     nuclear transparency, no formation-time physics. The partners are put on
     shell and never touched again.
   * **Triton remnant realism.** `t* -> n + d` / `t* -> p + (nn)` is a
     sequential two-body model with the AME2020 separation energies and an
     isotropic S-wave relative direction, and the unbound nn pair is split at
     its virtual-state pole `1/|a_nn|` = 10.44 MeV. plans/05 5.D asked for
     exactly this ("crude, flagged"); a Faddeev/AV18 three-body triton wave
     function with a correlated (p, n, n) momentum distribution is the
     replacement. Conservation does not depend on it -- the struck nucleon
     absorbs the whole difference -- so what the crudeness costs is the
     SHAPE of the partner spectra.
   * **No tensor structure in the triton breakup** (it is isotropic), where
     the deuteron's D wave is fully correlated with m_S.
3. **CLOSED (2026-08-30) -- coherent T2.** The gap named in the previous
   revision ("PythiaBridge v0 has no coherent-diffractive target at all";
   the inclusive fallback broke whole-record momentum by ~16% and charge on
   about half of events) is closed by the **gamma*-Pomeron tier**
   (`docs/PYTHIA_BRIDGE.md` §12, design `docs/open_items/code_designs.md`
   §1 option a'): `Pipeline::make_coherent` writes the T2 target itself --
   `Role::Pomeron`, `P_IP = P_ion - P_recoil`, pdg 990, status 3 -- and
   `PythiaBridge` hadronizes the gamma*-Pomeron system on a third PYTHIA
   instance (`Beams:idA = 990`, meson-like beam, antiquark remnant) through
   the same surrogate + frame map, with `W^2 -> M_X^2` and `zeta = beta`
   exactly. The old refusal (`hadronize_coherent`) is gone; the knob that
   replaced it is `PythiaBridgeOptions::coherent_t2 = {Pomeron, Off}` (CLI
   `--coherent-t2 pomeron|off`). Measured over 300 events (`test_t2.cpp`):
   300/300 hadronized, worst whole-nucleus 4p residual 4.8e-14 relative,
   charge exactly 0, |M_had - M_X|/M_X <= 2.3e-11, the intact recoil
   bit-identical, PYTHIA veto 0 at M_X >= 1.4 GeV (the M_X floor is 1.2,
   raised for exactly this tier), byte-identical HepMC3 across
   identical-seed runs. With `CoherentT2::Off` the bridge refuses every
   coherent event and the T0 record still closes exactly -- the
   `Role::HadronicX` pseudo-particle stays `Status::Final` and carries the
   whole system, so turning the tier off costs fidelity, never
   conservation. What remains open is physics, not bookkeeping: no
   exclusive-vector-meson channel below M_X = 1.2 GeV (design §1 option c,
   phase 2), and the Pomeron-PDF model choice (`PDF:PomSet`, default 6 = H1
   2006 Fit B LO) is the tier's largest systematic.
4. **CLOSED (2026-08-30) -- the "no surrogate" tail.** It was ~12 % for
   7Li-alpha and ~5-7 % for 6Li with the whole-cluster hook, because a whole
   off-shell cluster fed to the surrogate exceeded
   `PythiaBridgeOptions::headroom`. At T1 the surrogate is handed one
   nucleon, and the tail is **0 in 300 events on each of the three tagged
   channels**. `headroom` was not touched.
