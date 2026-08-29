# The PYTHIA 8 bridge (tier T2)

`lipolgen::PythiaBridge` (`include/lipolgen/pythia_bridge.hpp`,
`src/pythia/`) takes a T0/T1 event — beams, a scattered electron, a struck
nucleon or cluster — and appends the showered, hadronized final state of the
γ*–nucleon system as `Role::Hadron` particles.  It links stock **PYTHIA
8.317**; nothing in PYTHIA is patched.

```
Event (T0/T1)                    PythiaBridge                    Event (+T2)
 beams, e', p_N   ──►  surrogate e+N at (W², Q²)  ──► PYTHIA ──►  Role::Hadron
                        LHAup, frameType 5              │           Role::StruckNucleon
                                                        └─► frame map back
```

---

## 1. Why an `LHAup`, and why a surrogate

Two facts about PYTHIA 8.317 fix the whole design.

**(1) `Beams:allowMomentumSpread` runs the hard process at the *initial* √s.**
`ProcessLevel.cc:645` gates `PhaseSpace::newECM` on `doVarEcm`
(`Beams:allowVariableEnergy`), and `Pythia.cc:662` aborts if `doVarEcm` is
combined with a hard process.  With momentum spread the event is generated at
the initialisation √s and merely *boosted* into the actual beam directions.
For LHC beam spread (10⁻⁴) that is irrelevant; for a Fermi-smeared nucleon it
is not — at 100 GeV/u a 0.2 GeV Fermi momentum boosted by γ ≈ 107 spans
~81–124 GeV per nucleon, a ±20 % swing in √s.  `BeamSetup.cc:645` additionally
forces the beam on shell at the *PDG* mass, so a genuinely off-shell bound
nucleon is not representable either.  The only mechanism with per-event
four-momentum freedom is an in-memory `Pythia8::LHAup` with
`Beams:frameType = 5` (`docs/surveys/pythia8_survey.md` §4(v), §7).

**(2) Even through `LHAup`, PYTHIA forces the total final-state four-momentum
to the *initialisation* beam total.**  `BeamRemnants::setKinematics` sets the
remnant light-cone budget from `infoPtr->eCM()`
(`BeamRemnants.cc:662`: `wPosRem = eCM; wNegRem = eCM;`, and `:935`:
`wPosRem = eCM - (pSumOut.e() + pSumOut.pz())`), and `eCM` is fixed at `init()`
for a hard process.  Measured directly (50 events, LHAup with a target
20 GeV above nominal and 0.2 GeV of p_T): the total final-state four-vector
came out **identical in every event** and equal to the init beam total, not to
`k + p_N`.

Since √s is an *invariant*, no Lorentz transformation can reconcile a
per-event `(k + p_N)²` with a fixed `eCM²`.  So the bridge does not try.

### The surrogate

The hadronic final state of a DIS event is a property of the γ*–nucleon
subsystem: given `W²`, `Q²`, the struck flavour and its light-cone fraction,
the final state in the γ*N rest frame is fully specified.  The bridge
therefore

1. computes the *physical* invariants from the event: `q = k − k′`,
   `Q² = −q²`, `W² = (q + p_N)²`;
2. constructs a **surrogate** e + N event on PYTHIA's own fixed beams with the
   **same `W²` and `Q²`** — an ordinary collinear DIS event that PYTHIA is
   entirely happy with;
3. runs `pythia.next()`;
4. maps the resulting hadronic system onto the physical one with a **pure
   Lorentz transformation** (boost + rotation), fixed by aligning the γ*
   direction and the lepton plane in the two rest frames.

Because `W²` matches by construction, step 4 is a genuine Lorentz
transformation and distorts nothing: every hadron keeps its mass, the event
shape in the γ*N rest frame is untouched, and

    Σ hadrons + e′ = k + p_N

holds **exactly** (measured worst case 8.4 × 10⁻¹⁴ relative over 400 events).

#### Surrogate construction, in the PYTHIA lab frame

Beam A is the nucleon, `P_A = (E_A, 0, 0, p_A)` with `m_A` the PDG mass of
2212/2112; beam B is the electron, `k_B = (E_B, 0, 0, −E_B)`.  Solve for a
massless `k̃`:

    k_B · k̃ = Q²/2                       →  Ẽ′ + k̃_z = Q²/(2 E_B)
    P_A · k̃ = P_A · k_B − B̃              with B̃ = (W² − m_A² + Q²)/2

so that `q̃ = k_B − k̃` has `−q̃² = Q²` and `(q̃ + P_A)² = W²`:

    Ẽ′   = [ P_A·k_B − B̃ + p_A Q²/(2E_B) ] / (E_A + p_A)
    k̃_z  = Q²/(2E_B) − Ẽ′
    k̃_T  = sqrt(Ẽ′² − k̃_z²)                (azimuth 0; the map restores φ)

`E_A` and `E_B` are a **pure frame choice**: `W²`, `Q²`, the parton fraction
and the Bjorken `x` of the surrogate do not depend on them.  They only have to
be large enough that `k̃_T²  ≥ 0`.  `PythiaBridgeOptions::headroom`
(default **1.5**) sets `E_A = headroom × E_N,nominal`, which covers the ±20 %
Fermi tail with margin.  An event that still does not fit is counted in
`PythiaBridgeStats::n_no_surrogate` and `hadronize` returns false.

---

## 2. The ξ solution

Two fractions appear, and the header keeps them apart.

**Physical ξ** — the light-cone fraction of the *actual* (moving, possibly
off-shell) nucleon that puts the outgoing quark on its mass shell:

    (q + ξ p_N)² = m_q²
    ⇒ ξ² p_N² + 2 ξ (q·p_N) − (Q² + m_q²) = 0
    ⇒ ξ = (Q² + m_q²) / [ (q·p_N) + sqrt( (q·p_N)² + p_N² (Q² + m_q²) ) ]

This closed form (`lipolgen::dis_parton_fraction`) is chosen over the textbook
quadratic formula because it is numerically stable, it reduces **exactly** to
`ξ = x_Bj` for a massless on-shell target and a massless outgoing quark, and
for `p_N² < 0` (a strongly spacelike off-shell nucleon) it is the smaller of
the two positive roots.  It fails (returns false) when the discriminant is
negative or when `ξ > 1`.  This is the quantity reported by
`PythiaBridge::last_xi_physical()`.

**Surrogate ζ** — what PYTHIA actually receives.  The struck parton is a
massless collinear parton of beam A,

    p_in = (ζ P_A⁺ / 2) (1, 0, 0, 1),      P_A⁺ = E_A + p_A ,

exactly massless (so PYTHIA's `E = sqrt(m² + p²)` recomputation in
`ProcessContainer.cc` is a no-op) and exactly along +z (which
`LesHouches:matchInOut` requires).  ζ follows from the same mass-shell
condition:

    (q̃ + p_in)² = m_q²   ⇒   ζ = (Q² + m_q²) / (P_A⁺ q̃⁻)

and is handed over as `x1` in `setIdX`, which is then the true light-cone
fraction of beam A — PYTHIA's beam bookkeeping and beam remnant are exactly
consistent.  `PythiaBridge::last_xi_pythia()`.

ξ and ζ differ only through `p_N² ≠ m_A²` and the `m_A²q̃⁺` term dropped by
the light-cone definition; sub-percent at EIC energies.  ζ is the one used for
the PDF lookup, so the flavour choice and PYTHIA's backward evolution see the
same `x`.

**`m_q` is not always zero.**  `LesHouches:setQuarkMass = 1` (PYTHIA's
default) reassigns the table mass to a final-state **c or b** quark handed
over massless — and takes the recoil **from the scattered lepton**, which
would silently move the `e′` the core generator fixed.  Measured: this hit
8.8 % of events and shifted the hadronic `W` by up to 25 %.  The bridge
therefore puts c and b on their PYTHIA table mass in the ξ/ζ solve and in the
LHA record (`m_q = particleData.m0(|id|)` for |id| ∈ {4, 5}, zero for u, d,
s), so `setQuarkMass` finds nothing to fix.  After this the frame map is a
pure Lorentz transformation in every event: worst |λ − 1| = 1.2 × 10⁻¹³ over
400 events.

---

## 3. Flavour choice

The struck flavour is sampled with probability

    P(q) ∝ e_q² · x f_q(ζ, max(Q², q2_pdf_min))

using **PYTHIA's own beam PDF**, obtained from `pythia.getPDFPtr(idA)` once
per instance.  This is deliberately *not* routed through `sf.hpp`: `sf.hpp`
is a structure-function interface (F₁, F₂, g₁, b₁), not a PDF interface, and
using PYTHIA's own object guarantees that the flavour we pick and the PDF
PYTHIA's backward evolution starts from are the same grid at the same `x` and
scale.  The default is PYTHIA's NNPDF2.3 QCD+QED LO set.

Flavours offered: u, d, ū, d̄ always; s, s̄ with `include_strange` (default
on); c, c̄ with `include_charm` (default on); b, b̄ with `include_bottom`
(default **off**, since b is negligible in the EIC window and only slows the
sampler).

`q2_pdf_min` (default 1 GeV²) floors the PDF scale: the generator window
reaches `Q² = 0.7 GeV²` while the LO grid starts at 1 GeV².  Only the
*relative* flavour weights are affected.

**The neutron is not an isospin relabelling of the proton run.**
`PDF::resetValenceContent` (`PartonDistributions.cc`) sets `beamType = -1` for
2112, and `PDF::xf` then swaps u ↔ d and ū ↔ d̄ — so a neutron PDF comes for
free once the beam id is 2112.  More importantly the *beam remnant* depends
on the beam id: striking a d quark leaves a `uu` diquark in a proton and a
`ud` diquark in a neutron, with different charge, flavour and mass.  Beam ids
are `init()`-time settings, so **the bridge keeps two PYTHIA instances**, one
with `idA = 2212` and one with `idA = 2112`, and routes each event to the
right one.  (`with_neutron_instance = false` builds only the proton one.)
Both share one `RndmEngine` and both are initialised in ~0.34 s.

---

## 4. The LHA record

Per event, `LhaupDis::setEvent` writes the four-line LHEF process
(`src/pythia/lhaup_dis.hpp`):

| # | particle | LHA status | colour | momentum |
|---|---|---|---|---|
| 1 | struck quark | −1 | (101, 0) or (0, 101) | `ζ P_A⁺/2 · (1,0,0,1)` |
| 2 | electron | −1 | — | `k_B` |
| 3 | electron | +1 | — | `k̃` |
| 4 | quark | +1 | same tag as #1 | `q̃ + p_in`, mass `m_q` |

The incoming quark must come **first**: `ProcessContainer.cc` assigns the
first status −21 parton to beam A.  The colour tag runs from the incoming to
the outgoing quark, so the string is struck-quark ↔ beam remnant, which is
what makes `SpaceShower:dipoleRecoil` do the right thing.

`setIdX(id_q, 11, ζ, 1.0)` and `setPdf(..., √Q², xf/ζ, 1, true)` are set;
`LHAup` strategy 3 (unit weights, local process choice) is used, so PYTHIA's
`sigmaGen` is bookkeeping only and **not** a physics number — the cross
section is the core generator's (`Event::xsec_pb`).

`setEvent` stays armed after a call, because `Pythia::next()` may re-read the
same hard process when the parton level has to be retried.

**SPINUP.**  The incoming electron carries `Event::spin.lam_e` (±1, or the
beam particle's own `pol` when `lam_e == 0`); the incoming quark carries the
`Particle::pol` of the `Role::StruckNucleon`, i.e. the *nucleon* polarization
label, `9` when unknown.  Nothing in PYTHIA's DIS matrix element consumes it
(the survey §5: `SigmaProcess` is helicity-averaged); it is carried so the
record is honest and so τ decays and the weak shower see it if they ever fire.

---

## 5. The frame map

Let `H_p = Σ (PYTHIA final) − e′_PYTHIA` and `H_w = q + p_N`.  In each rest
frame build a right-handed triad

    ê₃ = γ* direction,   ê₁ = the incoming-lepton component ⊥ ê₃,   ê₂ = ê₃ × ê₁

(`ê₁` fixes the lepton plane, so whatever azimuthal correlation PYTHIA
produced about the photon axis is carried over faithfully), then for each
hadron: boost into the `H_p` rest frame, re-express the three-momentum in the
`H_w` triad, boost out into `H_w`.  Because the rotation is linear and each
rest-frame momentum sum is zero, the transformed sum is exactly `H_w`.

A common rescale λ solving `Σ sqrt(m_i² + λ²|p_i|²) = W_want` sits between the
two boosts to absorb any residual mass mismatch (safeguarded Newton).  With
the c/b mass fix and `SpaceShower:QEDshowerByQ = off` it is **1 to within
10⁻¹³** in every event measured; it is kept as a guard so that a settings
override cannot break conservation, and `PythiaBridge::last_rescale()` /
`PythiaBridgeStats::max_rescale_dev` report it.  A λ visibly different from 1
means PYTHIA moved the scattered electron and the event was repaired rather
than mapped — worth a look.

**The scattered electron** is identified as the bottom copy
(`Particle::iBotCopyId`) of the status-23 outgoing lepton of the hard process,
never by "most energetic final electron" (that is only a fallback that never
fires with the bridge's settings).  PYTHIA's copy is dropped; the Event's own
`Role::ScatteredElectron` is what stays, exactly as the core generator made
it.  The test checks that no `Role::Hadron` electron duplicates it.

**Spectators** (`Role::Spectator`, `Role::PartnerSpectator`) are never read and
never touched — no recoil correction, the BeAGLE light-nucleus rule
(`docs/surveys/beagle_survey.md`).

---

## 6. Choosing the target

| the event contains | what the bridge uses |
|---|---|
| `Role::StruckNucleon` | its four-vector and pdg, as given (may be off shell and moving) — **the only branch a `Pipeline` event ever takes** |
| `Role::StruckCluster` | **DEPRECATED** (2026-08-30, superseded by the T1 tier).  `p_cluster/A_c` on shell, flavour `Z_c : N_c`, no Fermi smearing, and **nothing carries the rest of the cluster**, so the whole record does not conserve: 0.186 relative and one charge unit wrong on half the events.  Warns once per bridge and counts into `PythiaBridgeStats::n_cluster_fallback`.  `set_nucleon_in_cluster(...)` still overrides it. |
| neither (inclusive T0) | a nucleon at rest in the ion rest frame, `p_N = P_ion/A`, on shell at the FREE `M_NUCLEON` (docs/CONVENTIONS.md).  Flavour `Z F2p : N F2n` by default (`NucleonChoice::ByStructureFunctions`), forced with `NucleonChoice::Proton/Neutron`, or replaced entirely with `set_nucleon_chooser(...)`. |

**Since the T1 tier, the cluster branch is unreachable from `Pipeline`.**
`PipelineConfig::tier` defaults to `Tier::T1`, which resolves the tagged
struck cluster into a struck NUCLEON plus its on-shell partner spectator(s)
(`breakup.hpp`) before the hadronizer hook runs, so `Role::StruckNucleon` is
always present and the whole record conserves with no hook at all
(`docs/T2_CHAIN.md` §1a).  A non-zero `n_cluster_fallback` on a pipeline run
means the run is at `Tier::T0`.

When the target was implicit the bridge **appends** the nucleon it used as a
`Status::Intermediate`, `Role::StruckNucleon` particle, so the record says
what was struck.

Any `Role::HadronicX` pseudo-particle is demoted to `Status::Intermediate`;
the hadrons that replace it are `Status::Final`.

---

## 7. Exact PYTHIA settings

Applied in this order, before any user `PythiaBridgeOptions::settings`
(which therefore win); `PythiaBridge::applied_settings()` returns the list
verbatim.

```
Beams:frameType = 5
Check:beams = on
SpaceShower:dipoleRecoil = on
SpaceShower:pTmaxMatch = 2
PDF:lepton = off
TimeShower:QEDshowerByL = off
SpaceShower:QEDshowerByQ = off
LesHouches:setLeptonMass = 0
Random:setSeed = on
Random:seed = <1 + seed mod 9e8>
Next:numberCount = 0
Next:numberShowInfo = 0
Next:numberShowProcess = 0
Next:numberShowEvent = 0
Print:quiet = on                    (verbosity 0 only)
```

- `SpaceShower:dipoleRecoil`, `SpaceShower:pTmaxMatch = 2`, `PDF:lepton = off`
  and `TimeShower:QEDshowerByL = off` are exactly the four DIS shower settings
  of `PolarizedLithiumSim/tools/pythia8/gen_dis_hfs.py` (= PYTHIA's
  `examples/main341.cc`).
- **`Check:beams = on`, not off.**  The survey suggested `Check:beams = off`
  for a nonstandard beam combination; that is backwards here.
  `BeamSetup::checkBeams` (`BeamSetup.cc:878-885`) accepts a lepton+hadron
  pair only when a `WeakBosonExchange:*` flag is on, **or `Check:beams` is
  on**, or `frameType == 4`.  We want neither an internal DIS process nor
  LHEF-from-file, so the flag has to be on (it is PYTHIA's default; the bridge
  sets it explicitly so a user override cannot silently break `init()`).
  Setting it *off* makes `init()` fail with
  "cannot handle this beam combination".
- **`SpaceShower:QEDshowerByQ = off`** — a QED ISR photon off the struck quark
  takes its recoil from the scattered lepton (the only other charge in the
  dipole), which would move the `e′` the core generator fixed.  Measured over
  500 events: 0.4 % of them, with the hadronic `W` shifted by up to 25 %.
  Lepton-side QED radiation is out of scope for v0.1
  (`DEVELOPMENT_PLAN.md` §2) and PYTHIA cannot do it consistently with dipole
  recoil anyway, so every QED emission that can reach the lepton leg is off.
  `TimeShower:QEDshowerByQ` stays **on**: photons radiated by the final-state
  quark belong to the hadronic system and were never observed to move the
  lepton.
- `LesHouches:setLeptonMass = 0` keeps the outgoing lepton exactly massless as
  handed over; the default (1) reassigns `m_e` and shuffles the difference
  onto the struck quark.
- `Beams:allowMomentumSpread` is **not** used, for the reason in §1.
- **No `PhaseSpace:*` setting is applied, and none would do anything.**  The
  hard process arrives through `LHAup`, so `PhaseSpace:Q2Min`,
  `PhaseSpace:pTHatMinDiverge` and `PhaseSpace:mHatMin` — including the two
  silent cuts that `gen_dis_hfs.py` had to discover and fix — are irrelevant
  to this tier.  The generator window (`Q² ≥ 0.7`, `y ∈ [0.004, 0.985]`,
  `W² ≥ 8`) is enforced by the core sampler, not by PYTHIA.
- Randomness: PYTHIA draws through a `Pythia8::RndmEngine` backed by the
  per-event `lipolgen::Rng`, so the same `(seed, run, bunch, event)` gives the
  identical event whatever the order or the thread count
  (`docs/CONVENTIONS.md`).  `Random:seed` only seeds the fallback stream used
  during `init()`.

---

## 8. HFS truth identities

`hfs_summary(ev)` sums over `Role::Hadron`, `Status::Final`:

    Σ (E − p_z),  Σ p_T (as a vector and its modulus),  Σ E, Σ p_z,
    Σ charge,  multiplicities (total / charged / neutral).

Two identities are pinned by the tests, both **including the beam remnant**:

1. `Σ (E − p_z)_hadrons = 2 E_e y + m_N²/(E_N + p_Nz)`
   (`hfs_sigma_empz_truth`).  This is the massive-target version of the
   massless `Σ = 2 E_e y` that the toy generator satisfies; the `m_N²/(E_N+p)`
   term is 4.4 MeV at 99.5 GeV/u (10.8 / 4.4 / 3.2 MeV at 40.8 / 99.5 / 137.5),
   which is exactly the offset `tools/pythia8/README.md` measured on the Python
   samples.  It is itself approximate at the level
   `Q² m_N²/(2 E_e (E_N + p_Nz)²)` — 1.1 × 10⁻⁵ GeV at `Q² = 10` — because
   `k′⁻` carries a residual target-mass term.  Measured worst relative
   deviation over 400 events: **8.7 × 10⁻⁵**, tested at rtol 10⁻³.
2. `|Σ p_T,hadrons| = p_T,e′` for a collinear target, and in general
   `Σ p_T,hadrons = p_N,T − p_T,e′`.  Measured worst deviation
   **8.3 × 10⁻¹³ GeV**, tested at 10⁻⁸ GeV (the Python project measures
   ≤ 5 × 10⁻¹¹).

`hfs_sigma_empz_exact(k, p_N, k′) = (k + p_N − k′)⁻` is the statement that is
exact for any target and is checked at rtol 10⁻¹⁰.

---

## 9. What is NOT modelled

- **Lepton QED radiation.**  `PDF:lepton = off`, `TimeShower:QEDshowerByL =
  off` and now `SpaceShower:QEDshowerByQ = off`: the lepton leg does not
  radiate and nothing recoils against it.  PYTHIA's own documentation
  (`SpacelikeShowers.xml:371-381`) states that QED radiation off the lepton is
  not implemented together with `dipoleRecoil`, so this is a PYTHIA limitation
  and not only a choice.  Radiative corrections are a separate work package
  (`DEVELOPMENT_PLAN.md` §2, plans/07 WP4).
- **Nuclear effects on the hadronic final state.**  No shadowing, no EMC
  modification of the *fragmentation*, no nuclear transparency, no FSI of the
  produced hadrons in the residual nucleus, no formation-time physics.  The
  nucleus enters only through (a) which nucleon is struck, (b) that nucleon's
  four-momentum.  Nuclear PDFs are hard-process only in PYTHIA and give an
  *average* nucleon (`PDFSelection.xml:333-345`), which is the opposite of
  what a tagged generator needs, so they are not used.
- **The remnant of an off-shell nucleon.**  PYTHIA builds the beam remnant
  from an *on-shell* PDG-mass nucleon of the surrogate.  The physical
  off-shellness (`p_N² ≠ m_N²`) enters the invariants `W²` and `Q²` exactly —
  those are matched — but the remnant's own internal composition, its
  primordial k_T and its diquark mass are those of a free nucleon.  For the
  `|p_N² − m_N²| ≲ 0.1 GeV²` of light-nucleus Fermi motion this is a
  sub-percent effect on the remnant's light-cone fraction, but it is an
  approximation, not an identity.
- **Spin in the hard process.**  `SigmaProcess` is helicity-averaged (survey
  §5); the polarized cross section lives entirely in the core generator's
  weight.  SPINUP is carried, not consumed.
- **The `NucleonInCluster` v0 has no Fermi smearing** — and is now
  DEPRECATED and unreachable from `Pipeline`.  What replaced it is not a
  better hook but a better record: `breakup.hpp` (tier T1) draws the internal
  relative momentum from the cluster's own wave function, emits the partner
  spectator(s) on shell and hands the bridge a real struck nucleon.  Fermi
  smearing therefore reaches the hard process the same way it does on the
  inclusive channel — through `Role::StruckNucleon`'s four-vector — and the
  hook is left only for a caller who builds a record by hand.
- **MPI / multiple scattering**: absent by construction with a lepton beam.
- **Charm/bottom as massive initiators**: the *incoming* parton is always
  massless and collinear (the standard collinear-factorisation statement); only
  the outgoing c/b is put on its mass shell.

---

## 10. Measured performance and cross-checks

Single core (this machine), 10 × 99.5 GeV/u ⁶Li, generator window:

| quantity | value |
|---|---|
| initialisation (both instances) | 0.34 s |
| hadronization only | **33 300 events/s** |
| hadronization + HepMC3 Asciiv3 write | **6 200 events/s** (16.5 MB / 5000 events) |
| PYTHIA vetoes / flavour retries | 0 in 5000 events |
| worst relative 4-momentum deviation | 8.4 × 10⁻¹⁴ |
| worst \|λ − 1\| | 1.2 × 10⁻¹³ |

Mean charged multiplicity of the hadronic final state, e + p, `Q² ∈ [4, 30]`,
`x ∈ [0.005, 0.10]`:

| generator | ⟨n_charged⟩ |
|---|---|
| stock PYTHIA `WeakBosonExchange:ff2ff(t:gmZ)`, `gen_dis_hfs.py` settings, 1153 events in the window (of 2000) | **8.09** at ⟨x⟩ = 0.0312, ⟨Q²⟩ = 8.88 |
| `PythiaBridge`, 2000 events at that exact (x, Q²) | **7.87** |
| ratio | **0.97** |

The residual 3 % is expected: the stock run averages over the whole window
while the bridge sits at its mean point, and ⟨n_ch⟩ grows like `ln W²`.  The
test asserts only |ratio − 1| < 20 %.

---

## 11. Files

| file | contents |
|---|---|
| `include/lipolgen/pythia_bridge.hpp` | public, PYTHIA-free (pimpl): `PythiaBridge`, `PythiaBridgeOptions`, `PythiaBridgeStats`, `HfsSummary`, `dis_scattered_electron`, `dis_parton_fraction`, `hfs_summary`, `hfs_sigma_empz_truth/exact` |
| `src/pythia/lhaup_dis.hpp` | internal: the `LHAup` subclass and the `RndmEngine` that routes PYTHIA's randomness through `lipolgen::Rng` |
| `src/pythia/pythia_bridge.cpp` | the surrogate, the flavour sampler, the frame map |
| `tests/test_pythia.cpp` | 7 doctest cases; the head-on convention, the ξ solve, 200 events each on p and n, a Fermi-moving off-shell target, determinism, the target-choice hooks, and the stock-PYTHIA multiplicity comparison |
| `examples/hadronize_example.cpp` | N synthetic events → HepMC3 + timing + bookkeeping |
