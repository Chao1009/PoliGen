# HepMC3 output convention

`lipolgen::HepMC3Writer` (`include/lipolgen/hepmc_writer.hpp`,
`src/hepmc/hepmc_writer.cpp`) writes `lipolgen::Event` records
(`include/lipolgen/event.hpp`) as HepMC3 Asciiv3 files. This document is the
**proposed convention for ion-spin states in HepMC3** flagged as open item
**#17** in `PolarizedLithiumSim/plans/04_open_questions.md` ("HepMC3
convention for ion spin states — none exists; plans/05 step 5.D defines named
attributes and proposes them upstream. *Engage:* ePIC MC/software group.
*Default:* our attribute schema."). No such convention exists upstream in
HepMC3 or in the EIC software stack, so LiPolGen defines one here and is
prepared to propose it to the ePIC MC group.

## Format

Only `HepMC3::WriterAscii` (Asciiv3) is implemented (`HepMC3Format::Asciiv3`,
the default and only working value of the `HepMC3Writer` constructor's format
argument; `HepMC3Format::HepMC2Ascii` is declared for API completeness only
and the constructor throws `std::runtime_error` if it is requested).

This is a hard requirement, not a style choice: DD4hep routes a `.hepmc` /
`.hepmc3` input file to its `HEPMC3FileReader`, which opens the file with
HepMC3's own `ReaderAscii`. Handing that a HepMC2 `IO_GenEvent` file produces
an immediate EOF ("Error when moving to event - EOF") rather than a
recognizable parse error — this was one of the harder-to-diagnose lessons of
`PolarizedLithiumSim/tools/fullsim/ion_gun_hepmc.py`. So: always Asciiv3.

Units are GeV / mm (`HepMC3::Units::GEV`, `HepMC3::Units::MM`), matching the
head-on frame convention of `docs/CONVENTIONS.md`.

## Event graph

- **Primary vertex.** A single `GenVertex` with both beam particles —
  electron and ion — incoming, status 4 (`HepMC3::Status` 4 = incoming beam
  particle, matching `lipolgen::Status::Beam`). The ion carries a 10-digit
  nuclear PDG code (`1000030060` = ⁶Li, `1000030070` = ⁷Li, per the standard
  `10LZZZAAAI` nuclear-code scheme; ⁴He spectators use `1000020040`).

  This vertex is always written normally — HepMC3 requires at least one
  vertex per event, and unlike a bare single-particle "gun" event with
  **nothing** incoming (where HepMC3's own `WriterAscii` silently omits the
  vertex line and a written `V` line with zero incoming particles fails
  `ReaderAscii`'s parser — see the `ion_gun_hepmc.py` docstring in
  `PolarizedLithiumSim` for the exact failure mode), our primary vertex
  always has two incoming beam particles, so it is emitted as an ordinary `V`
  line with no special-casing needed.

- **Everything else** is attached using `Particle::mother1`/`mother2`
  (indices into `Event::particles`) when set: a daughter is added as
  outgoing from the *end vertex* of its mother — reusing the mother's
  existing end vertex if it already has one (e.g. a particle whose mother is
  a beam particle lands back on the primary vertex, since the beam particle's
  end vertex *is* the primary vertex), or creating a new vertex with the
  mother(s) incoming otherwise. Particles with no mother set
  (`mother1 == -1`) are hung directly off the primary vertex. This lets a
  T0-only event (no explicit struck-cluster particle) still produce a
  physically sensible two-vertex graph: e.g. scattered electron, virtual
  photon and the spectator all parented by a beam particle land on the
  primary vertex, while the pseudo-particle X, parented by the virtual
  photon, gets its own downstream vertex.

  Caveat for producers: `mother2`, if given, must not name a particle that
  is already consumed (has an end vertex) at a *different* vertex than
  `mother1`'s end vertex — most concretely, don't list a raw beam-particle
  index as a second mother once something else has already consumed that
  beam at the primary vertex. Route through an explicit intermediate
  particle (struck cluster/nucleon) instead, the way `StruckCluster` /
  `StruckNucleon` roles are intended to be used in T1/T2.

- **Status mapping.** `Particle::status` is written through as the HepMC3
  status code (`Status::Beam = 4`, `Final = 1`, `Decayed = 2`,
  `Intermediate = 3`), **except** `Role::HadronicX`, which is always forced
  to status 3 (documentation), regardless of the stored `Status` value. This
  guarantees X is never written as a final-state particle — in particular,
  when the event also carries real T2 hadrons (`Role::Hadron`, status 1) as
  its actual final state, X still appears (for provenance/bookkeeping) but
  is never double-counted as final. The virtual photon (`Role::VirtualPhoton`)
  and the coherent channel's Pomeron (`Role::Pomeron`, PDG 990,
  `P_IP = P_ion − P_recoil`, spacelike, so its `mass` follows the photon's
  convention and is written negative, `−√|t|`) are likewise expected to
  arrive already marked `Status::Intermediate` by the
  producer; documentation-only status-3 particles like them and X don't
  contribute to `Event::total_final()` and shouldn't be summed a second time
  when checking overall 4-momentum conservation on the HepMC3 side — see the
  "conservation" note below.

- **`Particle::mass`** is written via `GenParticle::set_generated_mass`, so
  it round-trips independently of any small numerical mismatch between the
  stored 4-vector and the physical on-shell mass. One exception: a massless
  electron (`|pdg| == 11`, `p.mass == 0.0`, the core's standard DIS
  kinematics) is written with `generated_mass = 0.51099895e-3` GeV (the PDG
  m_e) instead of 0, so downstream Geant4/DD4hep sees a consistent record
  and doesn't nudge E by O(10 ppm) to enforce E² − p² ≥ m_e².

- **Role → status/PDG map**, complete, with the tier that writes each row:

  | role | status | PDG | tier |
  |---|---|---|---|
  | `BeamElectron`, `BeamIon` | 4 | 11 / 10-digit ion | T0 |
  | `ScatteredElectron` | 1 | 11 | T0 |
  | `VirtualPhoton` | 3 | 22 | T0 |
  | `Spectator` | 1 | 10-digit ion (`1000020040` for the α), or 2212/2112 for a nucleon spectator | T0 |
  | **`PartnerSpectator`** | **1** | 2212 / 2112 for a nucleon, 10-digit ion for a bound remnant (`1000010020` for the t\* → n + d deuteron) | **T1** |
  | `StruckCluster` | 3 | 10-digit ion of the cluster (2212/2112 when the "cluster" is one nucleon) | T0 |
  | `StruckNucleon` | 3 | 2212 / 2112 | T0 (inclusive) / **T1** (tagged) |
  | `IntactRecoil` | 1 | 10-digit ion | T0 (coherent) |
  | `Pomeron` | 3 | 990 | T0 (coherent) |
  | `HadronicX` | forced 3 | 92 | T0 |
  | `Hadron` | 1 | PYTHIA's | T2 |

  **`PartnerSpectator` is a genuine final-state fragment, status 1**, exactly
  like `Spectator`: it is one of the nucleons (or the bound d / the two
  neutrons) the struck cluster broke into, it is on shell at its AME2020
  mass, and it carries the cluster's 10-digit ion code when it is a nucleus.
  It differs from `Spectator` only in provenance — `Spectator` is the tagged
  cluster the far-forward detectors measure and the one `route_of` /
  `rp_tagged` classify, `PartnerSpectator` is an interior fragment that no
  routing function looks at.  Its `mother1` is the `StruckCluster` it came
  from.  A consumer summing status 1 gets both, which is what the whole-record
  balance needs; a consumer reconstructing the *tag* must select
  `Role::Spectator`, not "every status-1 nucleus".

## 4-momentum conservation bookkeeping

`Event::total_final()` sums `Status::Final` (HepMC3 status 1) particles only.
When the hadronic final state is represented purely by the T0 pseudo-particle
X (status 3, not final), overall conservation must be checked as

```
sum(status == 1 particles) + X.momentum  ==  sum(beam momenta)
```

*not* `sum(status == 1) == sum(beams)` alone (X's momentum is real and not
otherwise represented), and *not* by also adding the virtual photon's
momentum on top of that (its momentum equals `beam_e - e'`, already implied
by the electron beam and the scattered electron both being present — adding
it again would double count). `tests/test_hepmc.cpp` checks exactly this
identity to `1e-9`.

## Event attributes (spin / run labels)

All of the following are `GenEvent`-level attributes (id 0), named exactly
as listed, HepMC3 type in parentheses:

| name | type | source | notes |
|---|---|---|---|
| `spin_J` | `DoubleAttribute` | `Event::spin.j` | ion spin (1 or 1.5) |
| `spin_M` | `DoubleAttribute` | `Event::spin.m_ion` | ion projection on the quantization axis |
| `struck_cluster_m` | `DoubleAttribute` | `Event::spin.m_struck` | struck-cluster projection; NaN if inclusive |
| `lam_e` | `IntAttribute` | `Event::spin.lam_e` | electron helicity, ±1 or 0 |
| `P_e` | `DoubleAttribute` | `Event::spin.pe` | electron polarization magnitude |
| `P_z` | `DoubleAttribute` | `Event::spin.pz` | ion fill vector polarization |
| `P_zz` | `DoubleAttribute` | `Event::spin.pzz` | ion fill tensor polarization |
| `spin_axis_theta` | `DoubleAttribute` | `Event::spin.theta_s` | quantization axis polar angle, head-on frame |
| `spin_axis_phi` | `DoubleAttribute` | `Event::spin.phi_s` | quantization axis azimuth, head-on frame |
| `spin_category` | `StringAttribute` | `Event::spin.category` | `SpinCategory` name |
| `run` | `IntAttribute` | `Event::spin.run` | bookkeeping run id |
| `bunch` | `IntAttribute` | `Event::spin.bunch` | bookkeeping bunch id |
| `channel` | `StringAttribute` | `Event::channel` | enumerator name, e.g. `"TaggedLi6Alpha"` |
| `dis_x` | `DoubleAttribute` | `Event::kin.x` | |
| `dis_Q2` | `DoubleAttribute` | `Event::kin.q2` | |
| `dis_y` | `DoubleAttribute` | `Event::kin.y` | |
| `dis_phi` | `DoubleAttribute` | `Event::kin.phi` | |
| `spectator_k` | `DoubleAttribute` | `Event::kin.k` | spectator momentum, spin frame |
| `spectator_cos_theta` | `DoubleAttribute` | `Event::kin.cos_theta_k` | |
| `spectator_phi` | `DoubleAttribute` | `Event::kin.phi_k` | |
| `alpha_s` | `DoubleAttribute` | `Event::kin.alpha_s` | light-front spectator momentum fraction |
| `pt_s` | `DoubleAttribute` | `Event::kin.pt_s` | |
| `t` | `DoubleAttribute` | `Event::kin.t` | coherent channel |
| `x_pom` | `DoubleAttribute` | `Event::kin.x_pom` | coherent channel |

Cross section (pb) is carried the standard HepMC3 way, via a
`GenCrossSection` attribute (`GenEvent::set_cross_section`) set from
`Event::xsec_pb` / `Event::xsec_err_pb` on every event, not as a named
LiPolGen-specific attribute. LiPolGen always writes a single-entry
`cross_sections`/`cross_section_errors` vector (index 0), but
`GenCrossSection::from_string` **pads that vector out to the event's weight
count on read-back**, duplicating index 0 into the extra slots (see
`GenCrossSection.cc`) — so a reader should only ever trust index 0, not the
vector's size, as LiPolGen's actual cross section.

`Particle::pol` (PYTHIA/LHEF `SPINUP`-style helicity label) is written as a
per-particle `"pol"` `DoubleAttribute` only when it is not the "unknown"
sentinel `9.0`.

## Event weights and `GenRunInfo`

`GenRunInfo` is built once, from the first event written, and reused
(the same `shared_ptr<GenRunInfo>`) for every subsequent event so HepMC3
does not warn about mismatched run-info objects:

- **Tool info**: one `GenRunInfo::ToolInfo` entry, `{name, version}` from the
  `HepMC3Writer` constructor arguments (default `"LiPolGen"` / `"0.1.0"`).
- **Weight names**: `"nominal"` at index 0, plus one `"spin_weight_<k>"`
  (1-based `k`) per entry of `Event::spin_weights`, only if that first
  event's `spin_weights` is non-empty. Per-event weight *values* are always
  written explicitly (`GenEvent::weights()`, size `1 + spin_weights.size()`
  for that event), so an event whose `spin_weights` count differs from the
  one that established the run's registered names still round-trips its
  values correctly by index — it just isn't fully reachable by name through
  `GenRunInfo::weight_index`.

## Why not omit the primary vertex here

The one prior LiPolGen-adjacent HepMC3 writer in this codebase family,
`PolarizedLithiumSim/tools/fullsim/ion_gun_hepmc.py`, hand-rolls Asciiv3 text
for a single beam-less ion and has to omit the vertex line entirely (a vertex
with **no incoming particles** is dropped by HepMC3's own `WriterAscii`, and
writing it anyway desyncs the vertex count in the `E` line from what
`ReaderAscii` expects). `HepMC3Writer` never hits that case: the primary
vertex here always has exactly two incoming beam particles, so
`HepMC3::WriterAscii` is used as-is (no hand-rolled ASCII) and the vertex is
written like any other.
