# LiPolGen

**Doubly polarized e + ⁶Li / e + ⁷Li DIS event generator for the EIC.**
C++17 core (spin-density kernel, cluster–spectator sampler, coherent channel),
stock **PYTHIA 8.317** as the hadronizer (through `LHAup`, no fork), **HepMC3**
output with an ion-spin attribute convention, and a **pybind11** Python module
that is a drop-in, ~100× faster replacement for the numpy generator
`PolarizedLithiumSim/evgen/polligen` it was ported from.

License: **GPL-3.0-or-later** (`LICENSE`). Every source file under
`include/`, `src/`, `tests/`, `python/lipolgen/`, `python/bindings.cpp` and
`validation/` carries an `SPDX-License-Identifier: GPL-3.0-or-later` line —
**101 of 101 files** at the last run of `validation/check_spdx_headers.py`,
which gates it. `AUTHORS` names the copyright holder(s) and `CITATION.cff`
says how to cite this software; **both currently carry the literal
placeholder `<AUTHOR NAME — to be filled by the author>` in their one name
slot** — no name was guessed, `python/tests/test_release_metadata.py` keeps
the two files in sync, and filling them in is an open author decision
(`docs/open_items/run_2026-09-03/STATUS.md` row 13). No copyright line was
added to any source file: the identifier carries no holder.

## What it generates

| channel | final state | physics |
|---|---|---|
| `inclusive` | e′ + struck nucleon (p/n ∝ Z F₂ᵖ : N F₂ⁿ) + X | HJM spin-1 / spin-3/2 master formula: F₁, F₂, R, g₁, g₂^WW, b₁, b₂, Δ (cos 2φ gluonometry); E143 finite-γ A∥ (default); exact finite-γ tensor kernel (Cosyn 2025, `tensor_gamma`, off by default) |
| `tagged-6Li-alpha`, `tagged-7Li-alpha`, `tagged-d-p` | + spectator fragment (α / α / p) with the **spin ⊗ cluster-wave-function** correlation; T1 breakup of the struck cluster (d → N + N, t → N + d/nn, or the three-channel **Ciofi–Simula** spectral function via `--triton-sf ciofi-simula`) | two-cluster light-front impulse approximation; radial forms: Hulthén S/P/D (default) or **ANL VMC AV18** α+d, α+t tables (`--cluster-wave vmc`); optional **Glauber spectator-FSI weight** (`--fsi`, a weight on `Event::weight`, never a momentum shift) |
| `coherent` | e′ + intact ⁶Li recoil + X | scenario f_coh(x)·e^{−B|t|} with deformation + gluon-transversity cos 2φ; x_P sampled with a physical M_X ≥ 1.2 GeV; T2 hadronizes γ*–Pomeron on a PYTHIA Pomeron beam (id 990, `--coherent-t2`) |

Every event carries its spin labels (λ_e, M, m_S, axis, P_z, P_zz, category,
run, bunch); run plans mirror `polligen.bookkeeping` (helicity flip, tensor
thirds, transverse tensor, tensor flip, relative-luminosity offsets, polarimetry).
Tier **T2** adds the PYTHIA hadronic final state on the *actual* Fermi-smeared,
off-shell struck nucleon (a (W², Q²)-matched surrogate γ*N event Lorentz-mapped
onto the physical target — conservation to 1e-13).

## What ships behind a selector

Almost every physics ingredient beyond the toy backends is **opt-in**: the
default run is the one every reference JSON is pinned against at rtol 1e-12,
and each selector below is off unless asked for. **The cost of turning one on
is measured, and no number here stands without the window it was measured in.**

| selector | default | what the other value adds | measured, with its window |
|---|---|---|---|
| `--cluster-wave {hulthen,vmc}` | `hulthen`, bit for bit | ANL VMC AV18 α+d / α+t tables, and the AV18 deuteron everywhere the run reads one | ⁶Li α-tag fraction **0.0249 → 0.0348** at 10 × 99.5 on the YR high-acceptance envelope; on the α-tag channel the opt-in path moves every polarized observable by **−2.027 %**; refused on `inclusive`/`coherent`, where it is never read |
| `--triton-sf {hulthen,ciofi-simula}` | `hulthen`, bit for bit | the Ciofi–Simula three-channel A = 3 spectral function | S₀ = **0.6525, untuned** |
| `--fsi {off,glauber-cluster,glauber-nucleon}` | `off` | a Glauber survival **weight** on `Event::weight` — never a momentum shift | the two variants are **not** one: **99.50 %** of events differ by more than **1 %** (\|w_nucleon/w_cluster − 1\| > 0.01), ratio to **68.5**, on `tagged-6Li-alpha`, 20 000 events, seed 1234, `tensor-thirds`, **σ_XN = 40 mb** — one end of the mandatory 20–40 mb band |
| `--rc {off,tensor-band}` | `off`, byte-identical | the tensor RC **band** `rc_tensor_lo/hi` plus the radiative tails, on `Event::rc_weights` and never on `Event::weight` | the band is the answer; a single edge is not. `--rc` is accepted and does nothing on `coherent`, and a non-default rc sub-knob there is refused |
| `--rc-tail-model {t-peak,t-peak+ll}` | `t-peak`, bit for bit | POLRAD's t-peak **plus** the leading-log s-/p-peaks | the two edges agree to **+0.61 % event-weighted** in the Q² ≥ 20 GeV², y ≤ 0.9 window and disagree **per cell**: 331 of 1356 accepted cells (24.4 %) by > 1 %, worst ×6444 at x = 0.7943, y = 0.0088 |
| `--rc-c0-shape {ho,vmc-ft}` | `ho`, bit for bit | the j₀ transform of the committed ANL VMC ⁶Li point-proton density at the same ⟨r²⟩ | the ⁶Li C0 shape is a **band**, and both the sign and the magnitude are band edges: σ^el_T/σ^el_U at x = 0.10, Q² = 5 is **+9.357e−04 on `ho`** and **−2.442e−04 on `vmc-ft`** — a factor **68** resp. **≈ 260** below POLRAD's own deuteron value +0.064, with the SIGN flipping between the edges, so the published sign is withdrawn |
| `--b1-model {miller,cdks,li6-convolution}` | `miller`, bit for bit | the four-term α–d convolution for b₁(⁶Li) | opt-in and **band-mandatory** (`--b1-band-scale` 0/1/2); the ±100 % band comes from Q(⁶Li) vs Q_d, not from the A = 2 gate, so it stays after the gate passes |
| `--b1-unpol {toy,mstw,ct18nlo}` | `toy` | MSTW2008 LO — CDKS's own PDF — read from the grid PYTHIA already ships | this is the selector that **emits the gate's passing configuration**: G3b **0.843243** at CDKS Eq. (21), inside [0.5, 2], against **0.440** on the shipped `toy`, outside it. Needs the optional PYTHIA tier and is refused, never silently downgraded, without it |
| `--unpol-sf {toy,mstw,ct18nlo}` | `toy`, bit for bit | one unpolarised structure-function backend for **every** kernel the run builds, the tagged struck-cluster kernel included | cost of having been on the toy, ⁶Li at config 1: accepted σ **×0.7985** (ct18nlo) / **×0.7934** (mstw), run-level A_zz ×1.2524 / ×1.2604 |
| `--pol-sf {toy,nnpdfpol}` | `toy`, bit for bit | NNPDFpol1.1 g₁ on the inclusive and tagged kernels | it is read **only** where the fill carries λ_e·P_e ≠ 0 — under `tensor-thirds`, this CLI's default plan, on **no** channel — and it is *labelled, not credited*, wherever it did not run. The shipped `ToyG1`'s g₁ⁿ has the **wrong sign over roughly 0.25 < x < 0.6** |
| `--coherent-t2`, `--pom-set` | `pomeron`, set 6 — both read only under `--hadronize`, itself off | the γ*–Pomeron hadronic final state on a PYTHIA Pomeron beam (id 990) | over all 15 sets × 20 000 events the T0 columns are **bit-identical across the fourteen sets that still run** (set 11 is refused), so the DPDF band on M_X, \|t\|, x_P and σ is identically zero; the systematic is hadronic — ⟨n_charged⟩ **−2.6 % / +9.9 %** over the twelve DPDF fits, kaons ×2.8 |
| `--coherent-t-max` | 0.2 GeV² | a different coherent \|t\| ceiling | the ceiling's reason is the **anchor range** (\|t\| ≤ 0.30, knob-independent); positivity is the contingent second reason — its edge is 0.245 at the shipped `eps_b0` and 2.80 at the measured quadrupole. Truncation redistributes exp(−B t_max) = **4.5e−5** of the rate |

**Every knob's read / not-read / refused status is one table** —
`Pipeline::knob_provenance`, **66 rows on a default run** — written into
`meta["knob_provenance"]` and printed as the banner's `KNOB PROVENANCE` block.
A knob that did not run is never recorded as if it had: the value is replaced
by a label (`not read on channel coherent-6Li`, `not read by plan
tensor-thirds`), and a value that would name *a variation of a piece that did
not run* is refused outright rather than recorded. The three far-forward
routing knobs are the one place where the label carries the value too — on a
channel that writes a fragment the classifier **is** consulted and only the
sensitivity is absent, so the label reads `not read on this run at 10x100
high-acceptance (…)` and the envelope name survives into `meta["optics"]` and
the banner header.

Event counts: `--events N` for a fixed count, or `--lumi L` **alone** for
luminosity mode (the two are exclusive, and giving both — on the command line
or through `--config-file` — is refused rather than silently resolved).

## Status (2026-09-05)

- 401 doctest cases / 17 240 286 assertions and 926 pytest cases pass, with
  1 doctest case and 112 pytest cases skipped. `pytest -rs` prints all 112
  with their reason — 28 lines, one per refused base configuration, each
  naming the cell and quoting `PipelineConfig::validate()`'s refusal
  (`python/tests/test_knob_provenance.py`, the knob-provenance matrix). The
  doctest skip is RC test T9, the `PolradFull` Rosenbluth limit, a code path
  v0 does not implement; its name and its reason are here and at
  `tests/test_rc.cpp:1609` and nowhere else, because doctest 2.4.11 reports a
  skipped case as a **count only** — a run filtered to that one case prints
  `402 skipped` and no name at all.
- A case that is nothing but an optional tier or an optional data tree is
  skipped at REGISTRATION time (`doctest::skip(...)`), so with an absent MSTW
  grid or `data/vmc` it is tallied SKIPPED rather than passed having measured
  nothing. Three rows are the stated exception, because each sits in a case
  whose shipped-default rows must still run: T16's `--b1-unpol mstw` reach row
  (`tests/test_b1_nuclear.cpp:2034`), T17's verdict-configuration row
  (`tests/test_b1_nuclear.cpp:2147`) and T1r's companion report
  (`tests/test_b1_nuclear.cpp:892`). Each prints a `SKIPPED (MSTW2008 LO
  unavailable)` message inside a case that is still tallied as **passed** —
  the open decision is registry row 25
  (`docs/open_items/run_2026-09-03/AUTHOR_DECISIONS.md` §B26).
- The strict citation gate reports **0 broken** over three documents: 1205 refs
  / 97 ranges / 7 external / 6 allow-listed on `PHYSICS_CHANNELS.md`, 19 refs /
  115 ranges on `theory/SPIN32_FINITE_GAMMA.md`, 8 external on
  `PYTHIA_BRIDGE.md`; `validation/check_spdx_headers.py` 101/101. All measured
  2026-09-05 at the end of the open-items run — the counts move with every
  phase, so read them off the gates and not off this line.
- Kernel, ρ-moments, tagged densities, spectator boosts, coherent scenario and
  bookkeeping agree with `polligen` (run 16, tensor sign to the literature
  convention) at **rtol 1e-12** against `validation/reference/*.json`.
- External anchors reproduced: Cosyn Eq. 27, Cosyn–Weiss deuteron TABLE II
  (+1/−2), Cosyn 2025 Table 1 finite-γ rows, ⁷Li ⟨P₂⟩ = −T/5, P_p = 0.866.
- ePIC chain gate passed: HepMC3 → `npsim` (direct, and via `abconv -p ip6_hiacc_100x10`).
- Throughput, single core: T0 inclusive 2.8 M ev/s, tagged 1.0 M, coherent
  2.3 M; +PYTHIA 33–44 k ev/s; HepMC3 writing ~7 k ev/s; Python columnar ~5×10⁵ ev/s.
- Independent adversarial review (`docs/code_review_2026-08-29.md`): all findings fixed.

Physics channels, models, references and implementing symbols: `docs/PHYSICS_CHANNELS.md`.

Open items, with explored solutions and prototypes: `docs/OPEN_ITEMS_SOLUTIONS.md`.
The 2026-09-03 open-items run that last revised all of it keeps its own
records in `docs/open_items/run_2026-09-03/`: **`SUMMARY.md`** (one page —
what the run measured and what it did **not** establish), `STATUS.md` (the
phase board and the ONE decision registry — cite a decision by its row
number), `AUTHOR_DECISIONS.md` (all 25 decisions stated in full with their
options and measured cost) and the per-phase number tables.
Rows 5–7, 9, 10, 12 and 14 now **ship code** (coherent T2 via a PYTHIA Pomeron
beam, `--coherent-t2`; the Ciofi–Simula triton spectral function,
`--triton-sf`; the Glauber spectator-FSI weight, `--fsi`; the tensor-sector
radiative-correction band and radiative tails, `--rc tensor-band`; the
four-term α–d convolution for b₁(⁶Li), `--b1-model li6-convolution` — opt-in,
band-mandatory, and since 2026-09-03 through its A = 2 validation gate, see
below; packaging, `pip install -e .`; one structure-function backend for
**every** kernel, `--unpol-sf toy|mstw|ct18nlo` / `--pol-sf toy|nnpdfpol` —
row 14, and the toy costs ×0.80 on σ(⁶Li) while the toy g₁ⁿ has the wrong sign
over 0.25 < x < 0.6). Rows 8 and 15 the board calls **deferred with design**:
row 8 is the spin-3/2 finite-γ theory note, which changes no behaviour
(`docs/theory/SPIN32_FINITE_GAMMA.md`); row 15 is the ⁷Li inclusive rank-2
sector, **exactly zero**, saying so since phase D in the run banner and in
`meta["rank2_input"]`, with the α–t b₁ that would fill it deferred on
decision row 20. Each item's state — closed, answered as a band, opt-in
shipped, deferred with design, or an author decision — is the board at the top
of that file, with the deciding number and its window.
Row 11's coherent-⁶Li amplitude on-ramp gained the eSTARlight unpolarized
baseline (`docs/open_items/run_2026-09-02/estarlight_li6.md`) and the polarized
α+d configuration sampler (`cluster_config.hpp`, console script
`lipolgen-configs`), and 2026-09-04 an answer: the ⁶Li tensor a₂ is
**marginally measurable** in coherent J/ψ at one EIC year — S = **2.63 σ**
at the band's low edge, the top a span **2.84 … 3.29 σ** whose crossing of
3 σ is *not* established, and 3 σ at **8.3 … 13.0 fb⁻¹/u** — from a
closed-form map, not a dipole-model amplitude, with no detection efficiency
below Q² = 0.1 in any source this tree has seen. Still open: the exact elastic radiative tail (Mo–Tsai,
beyond POLRAD's transcribed t-peak), the tensor cos 2φ dipole-model run with
the Mäntysaari group, b₁(⁶Li) beyond the α–d picture, and the polarized
quasi-elastic tail. **The A = 2 convolution gate now passes — for one
unpolarized nucleon input, not for the shipped default.** With MSTW2008 LO —
CDKS's own nucleon PDF, checklist item 4 — at CDKS Eq. (21)'s δ-function the
peak ratio is **0.843** against the digitized Fig. 4, inside the factor-2
window; on the **default** toy F₂ it is **0.440**, *outside* it. So quote
numbers made with `--b1-unpol mstw`, the selector that emits the passing
configuration (it needs the optional PYTHIA tier and is refused, never
silently downgraded, without it); the default `toy` backend's numbers are not
covered by the lift, and the two are not a rescaling of each other — mstw/toy
on b₁ is 1.848 / 1.276 / 0.817 at x = 0.10 / 0.30 / 0.50. With the real
CD-Bonn wave function as well the ratio is **1.000338**, i.e. a residual below
the error of digitizing a published figure — which is not the same claim as
three-digit agreement with CDKS, and it is specific to CD-Bonn *and* MSTW
together. The ⁶Li publication ban that failure imposed
is lifted; the mandatory ±100 % band is not, because it comes from Q(⁶Li) vs
Q_d and not from the gate. Of the three smaller decisions left open by the
close-out of the run before this one (`a94fd6e`, 2026-09-03), one was
**applied and awaits confirmation or reversal** — `CdksB1` stopped halving the
CDKS column, which is per nucleon already (decision registry row 1) — and two
are recorded as author decisions and still open: whether
`Li6ConvolutionOptions` should default to `r1998` like the A = 2 gate rather
than to `r_sigma_lt` (row 3), and Miller's own b₁ normalisation, which his
paper is self-inconsistent about by exactly the factor in question (row 2).
The registry is `STATUS.md`'s decision table; each row is stated in full, with
its options and their measured cost, in `AUTHOR_DECISIONS.md` (§B12, §B7,
§B2).

## Build

Dependencies (built from source into `../deps/install` on the development
machine; any prefix works): HepMC3 ≥ 3.2, LHAPDF ≥ 6.5 with `CT18NLO`,
`NNPDFpol11_100`, `EPPS21nlo_CT18Anlo_Li6`, PYTHIA ≥ 8.310 (with `--with-python`
if you want `import pythia8` too), pybind11 ≥ 2.10, Python ≥ 3.10 with numpy.

```bash
source env.sh                          # sets PATH/LD_LIBRARY_PATH/PYTHIA8DATA/LHAPDF_DATA_PATH/PYTHONPATH
cmake -S . -B build -DLIPOLGEN_DEPS_PREFIX=/path/to/prefix
cmake --build build -j8
./build/lipolgen_tests                 # C++ suite
(cd python && python -m pytest tests)  # Python suite
```

CMake options: `LIPOLGEN_WITH_LHAPDF`, `LIPOLGEN_WITH_HEPMC3`,
`LIPOLGEN_WITH_PYTHIA`, `LIPOLGEN_WITH_PYTHON`, `LIPOLGEN_BUILD_TESTS` (all ON).
Data tables (`data/vmc`) are found via `LIPOLGEN_DATA_DIR` (default: the source tree).

Alternatively, `pip install` the Python package on its own (scikit-build-core;
no `env.sh`, no `build/` needed):

```bash
LIPOLGEN_DEPS_PREFIX=/path/to/deps/install pip install -e .
```

`PYTHIA8DATA`/`LHAPDF_DATA_PATH` still need exporting at run time (see
`docs/USAGE.md`); the wheel's RPATH points at this machine's deps prefix, so
it is not relocatable as-is. `auditwheel repair` makes it mostly so —
measured 2026-09-05 with auditwheel 6.8.2 + patchelf 0.19.1: it vendors
`libpythia8`, `libLHAPDF` and `libHepMC3` in, **1 779 767 B → 7 232 971 B
(×4.06)**, tag **`manylinux_2_35_x86_64`** (glibc 2.35 — a wider tag needs a
build inside an older manylinux image). With the whole deps prefix replaced
by an empty tmpfs, `import lipolgen` and the pure-C++ generator **work**
(bit-identical output) while `LhapdfSF` and `--hadronize` **fail**, because
PYTHIA's `xmldoc`/`pdfdata` and LHAPDF's set store are not in the wheel and
their compiled-in defaults are the build machine's absolute paths. Vendoring
those three libraries also makes the wheel a combined work under
GPL-3.0-or-later. Nothing has been published; commands and measurements:
**`docs/PACKAGING.md`**, decision row 15.

`.github/workflows/ci.yml` (4 jobs) builds that whole dependency stack into a
cache and runs both suites and both gates against it (plus a PYTHIA-tier-off
build and the wheel). It has been exercised command by command on one Ubuntu
22.04 machine — deps 337 s at −j8, job 2 reproducing this machine's own
**400 cases / 17 240 262 assertions** bit for bit at the commit it was
written against, cache 431 MiB on disk / 113.5 MiB compressed — and it has
**never been run by GitHub**: nothing here parsed the file as GitHub parses
it, so cache round-trip, the runner image's apt set, real 4-core timings and
workflow-expression validity are unverified, and the file's own header says
so. It is committed and declares `on: push: branches: [master]`, so pushing
the run activates it — that coupling is decision row 14.

## Use

Command line (`./build/lipolgen-run` or `python -m lipolgen`):

```bash
lipolgen-run --isotope 6Li --config 1 --channel tagged-alpha --plan tensor-thirds \
             --events 100000 --seed 1 --optics tagging --cluster-wave vmc \
             --hepmc out.hepmc --npz out.npz --hadronize
```

Python:

```python
import lipolgen as lg
ev = lg.run(channel="tagged-6Li-alpha", isotope="6Li", config=1,
            plan="tensor-thirds", events=200_000, seed=1)   # columnar numpy dict
ev["k"], ev["cos_theta_k"], ev["R"], ev["m_ion"], ev["category"]
lg.export.write_hfs_npz(events, "hfs.npz")                   # polligen HFSSample format
```

C++: see `docs/USAGE.md` (`Pipeline`, `PipelineConfig`, `PythiaBridge`,
`HepMC3Writer`) and `examples/`.

## Layout

```
include/lipolgen/   spin sf xsec asymmetries beams bookkeeping sampler generator
                    spectator cluster tagged coherent breakup pipeline event rng
                    hepmc_writer pythia_bridge lhapdf_sf
                    rc b1_nuclear cluster_config fsi triton_sf mstw_sf
src/core src/lhapdf src/hepmc src/pythia
python/             bindings.cpp, lipolgen/{__init__,export,cli,configs}.py, tests/
tests/              doctest suites (one file per module) + reference-table comparison
validation/         dump_polligen_reference.py, reference/*.json, VMC reconciliation,
                    the two gates (check_physics_channels_links.py +
                    physics_channels_ranges.json, check_spdx_headers.py),
                    o5_a2_reach.py, o3_alpha_correlation_bound.py
data/vmc/           ANL VMC α+d / α+t overlaps and momentum distributions (provenance in README)
docs/               DEVELOPMENT_PLAN, CONVENTIONS, USAGE, HEPMC3_CONVENTION, PYTHIA_BRIDGE,
                    T2_CHAIN, PACKAGING, code review, OPEN_ITEMS_SOLUTIONS, PHYSICS_CHANNELS,
                    theory/ (SPIN32_FINITE_GAMMA), surveys/, open_items/ (design notes,
                    measured-number tables and gates per run)
LICENSE AUTHORS CITATION.cff   GPL-3.0-or-later, the copyright holder(s) and how to cite
                    (AUTHORS and CITATION.cff carry a name placeholder -- decision row 13)
.github/workflows/  ci.yml + scripts/build_deps.sh -- exercised locally, never run by GitHub
```

## Conventions (see `docs/CONVENTIONS.md`)

Head-on frame (ion +z, electron −z); populations ordered m = +J … −J;
`TENSOR_LL_SIGN = −1` (A_zz = −(2/3) b₁/F₁, Cosyn Eq. 27 / HERMES — the
value four sources give, and an author decision still open to confirm,
row 22); ⁶Li effective polarization **0.811228** whole-nucleus in the
cluster picture (0.811/3 per nucleon) — **no inclusive ⁶Li polarization may
be quoted without the band 0.81 … 0.91**, which brackets the ab-initio
six-body VMC 0.848;
`EmcBaseline::Epps21`; γ-matched beam energies (⁶Li 40.8 / 99.5 / 137.5,
⁷Li 40.8 / 99.5 / 117.9 GeV/u); counter-based RNG keyed by (seed, run, bunch,
event) — identical output for any thread count.

## Citing

PYTHIA 8.3 (Bierlich et al., SciPost Phys. Codebases 8 (2022)), HepMC3
(Buckley et al., CPC 260 (2021) 107310), LHAPDF6, and the physics inputs listed
in `docs/CONVENTIONS.md` and `docs/open_items/physics_literature.md`.
