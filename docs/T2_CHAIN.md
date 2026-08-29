# The T2 chain: Pipeline -> PythiaBridge -> HepMC3 -> ePIC

One page: how to run it, what a T2 record contains, what `abconv` / `npsim`
expect of the file, and the two things that do not fit yet.

```bash
source env.sh && cmake -S . -B build -DLIPOLGEN_WITH_PYTHON=OFF \
  && cmake --build build -j8 && ./build/lipolgen_tests
./build/generate_full --channel inclusive  --events 200000
./build/generate_full --channel 6Li-alpha  --events 200000 --plan azz
./build/generate_full --channel 7Li-alpha  --events 200000 --plan apar
./build/generate_full --channel coherent   --events 200000 --plan azz
./build/generate_full --channel inclusive  --events 200000 --no-hadronize   # T0 only
```

`--isotope 6Li|7Li` (inclusive only -- the tagged channels imply their own,
coherent is 6Li only), `--config 0|1|2` (energy point), `--seed`, `--out
FILE` (`""` = no HepMC3), `--no-hadronize` (skip PYTHIA, print T0
bookkeeping only). Each run prints sigma per spin category, the far-forward
tag fraction, the PYTHIA veto/retry/no-surrogate counters, the observed
whole-record conservation residual, and events/s.

## 1. What a T2 record contains

One `lipolgen::Event` per collision. `particles[0]` = beam e-, `[1]` = beam
ion, then T0/T1 (Pipeline) and T2 (PythiaBridge) particles interleaved by
the order each tier appended them.

| role | status | tier | HepMC3 |
|---|---|---|---|
| `BeamElectron`, `BeamIon` | `Beam` (4) | T0 | incoming at the primary vertex |
| `ScatteredElectron` | `Final` (1) | T0 | outgoing, untouched by the bridge |
| `VirtualPhoton` | `Intermediate` (3) | T0 | documentation only (optional) |
| `Spectator`, `PartnerSpectator` | `Final` (1) | T1 (tagged) | outgoing; **never read or written by the bridge** -- bit-identical before/after `hadronize()` |
| `IntactRecoil` | `Final` (1) | T1 (coherent) | outgoing; same "never touched" guarantee |
| `HadronicX` (pdg 92) | `Final` -> **demoted to `Intermediate` (3) by `hadronize()`** | T0 | documentation only once real hadrons exist, so a status-1 sum is never double-counted |
| `StruckNucleon` | `Intermediate` (3) | T0 (inclusive, explicit) or **appended by `hadronize()`** (implicit target) | documentation: what nucleon 4-vector/pdg actually went into PYTHIA |
| `StruckCluster` | `Intermediate` (3) | T1 (tagged) | documentation; read but not modified by the bridge |
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
2. **Tagged channels: the v0 `NucleonInCluster` default does not conserve
   the whole record.** `docs/PYTHIA_BRIDGE.md` section 6 documents v0 as
   "no Fermi smearing"; measured here, it is worse than that -- feeding
   PYTHIA only `p_cluster/A_c` (one nucleon's worth) with no particle
   carrying the rest of the cluster's momentum, plus a stochastic Z_c:N_c
   flavour draw that ignores the cluster's own integer charge, breaks
   `spectator + hadrons + e'` vs `beam_e + beam_ion` by **~15-20% in
   momentum on every event and by one charge unit on about half of them**
   (300-event 6Li-alpha measurement). `generate_full`/`test_t2.cpp` install
   `set_nucleon_in_cluster` — a public extension point named for exactly
   this, not a source patch — to hand PYTHIA the cluster's own off-shell
   four-vector *whole* and its exact charge, which restores exact
   conservation (measured ~1e-13 relative, exact charge) at the cost of a
   ~5-15% "no surrogate" tail (the whole off-shell cluster sometimes will
   not fit the PYTHIA-side surrogate's `headroom`) and losing the v0
   default's isospin sampling (this hook always assigns the cluster's own
   net charge, e.g. always "proton" for the Z=1 embedded deuteron). The
   complete fix is upstream of this file: the cluster's non-struck
   nucleon(s) should be emitted as `Role::PartnerSpectator` (the role
   already exists in `event.hpp`; nothing fills it for the alpha-tag
   channels today), so the bridge can strike one real nucleon and let the
   partner carry the rest on shell, the same way the deuteron-control
   channel's own spectator already works.
3. **Coherent: PythiaBridge v0 has no coherent-diffractive target at all.**
   A coherent event carries neither `Role::StruckNucleon` nor
   `Role::StruckCluster` (the diffractive system X is not a struck
   nucleon), so `hadronize()` always takes the plain inclusive fallback --
   a nucleon at rest in the ion frame -- which has nothing to do with the
   true (small, largely transverse) momentum transfer `P_ion - P_recoil`.
   The call succeeds mechanically (no crash, no excess vetoes: 300/300 in
   the measurement here) and the bridge's own internal HFS identities hold
   to the numerical floor (they are exact relative to *whatever* target it
   used), but whole-record momentum is off by **~16% and charge is wrong on
   about half of events**. Unlike the tagged case, **no public hook lets a
   caller supply the coherent target's four-vector** (`StruckNucleon`
   resolution is the bridge's first-priority branch, but nothing in
   `Pipeline::make_coherent` writes one, and there is no
   `NucleonInCoherent`-shaped extension point). This is reported here
   rather than worked around: it needs either a new public hook on
   `PythiaBridge` or, more fundamentally, a diffractive-dissociation model
   (Pomeron-exchange PYTHIA machinery, not the DIS-surrogate this bridge
   implements) -- both out of this file's scope (`src/` is owned
   elsewhere). `generate_full --channel coherent` prints the residual
   plainly rather than hiding it.
4. **7Li-alpha's "no surrogate" tail is larger than 6Li's** (measured ~12%
   at 2000 events vs ~5-7% for 6Li, both with the whole-cluster hook): the
   triton cluster's off-shell four-vector, fed whole, more often exceeds
   `PythiaBridgeOptions::headroom` (default 1.5x the nominal per-nucleon
   energy) than the lighter deuteron's does. A larger `headroom` would
   trade this for a slower `init()`; not tuned here.
