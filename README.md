# LiPolGen

**Doubly polarized e + ⁶Li / e + ⁷Li DIS event generator for the EIC.**
C++17 core (spin-density kernel, cluster–spectator sampler, coherent channel),
stock **PYTHIA 8.317** as the hadronizer (through `LHAup`, no fork), **HepMC3**
output with an ion-spin attribute convention, and a **pybind11** Python module
that is a drop-in, ~100× faster replacement for the numpy generator
`PolarizedLithiumSim/evgen/polligen` it was ported from.

License: **GPL-3.0-or-later** (`LICENSE`). Every source file under
`include/`, `src/`, `tests/`, `python/lipolgen/`, `python/bindings.cpp` and
`validation/` carries an `SPDX-License-Identifier: GPL-3.0-or-later` line
(`validation/check_spdx_headers.py` gates it). See `AUTHORS` for the
copyright holder(s) and `CITATION.cff` for how to cite this software.

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

## Status (2026-09-03)

- 363 doctest cases / 17.2 M assertions and 208 pytest cases pass.
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
Rows 5–10 and 12 are now **implemented** (coherent T2 via a PYTHIA Pomeron beam,
`--coherent-t2`; the Ciofi–Simula triton spectral function, `--triton-sf`; the
Glauber spectator-FSI weight, `--fsi`; the spin-3/2 finite-γ theory note,
`docs/theory/SPIN32_FINITE_GAMMA.md`; the tensor-sector radiative-correction
band and radiative tails, `--rc tensor-band`; the four-term α–d convolution for
b₁(⁶Li), `--b1-model li6-convolution` — opt-in, band-mandatory, and since
2026-09-03 through its A = 2 validation gate, see below; and packaging,
`pip install -e .`).
Row 11's coherent-⁶Li amplitude on-ramp gained the eSTARlight unpolarized
baseline (`docs/open_items/run_2026-09-02/estarlight_li6.md`) and the polarized
α+d configuration sampler (`cluster_config.hpp`, console script
`lipolgen-configs`). Still open: the exact elastic radiative tail (Mo–Tsai,
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
2026-09-03 close-out, one is now taken — `CdksB1` stopped halving the CDKS
column, which is per nucleon already — and two are recorded as author decisions
in `docs/OPEN_ITEMS_SOLUTIONS.md` §10: whether `Li6ConvolutionOptions` should
default to `r1998` like the A = 2 gate rather than to `r_sigma_lt`, and
Miller's own b₁ normalisation, which his paper is self-inconsistent about by
exactly the factor in question.

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
it is not relocatable as-is. `auditwheel repair` makes it so — measured
2026-09-05: it vendors HepMC3/LHAPDF/PYTHIA8 in, 1.70 MiB → 6.90 MiB, tag
`manylinux_2_35_x86_64`, and the result imports on a machine with no deps
prefix at all — but it pulls in GPL-3 redistribution terms for the combined
work, and the two libraries' **data** trees still are not in the wheel.
Commands and measurements: **`docs/PACKAGING.md`**.

`.github/workflows/ci.yml` builds that whole dependency stack into a cache and
runs both suites and both gates against it (plus a PYTHIA-tier-off build and
the wheel). It has been exercised command-by-command on one Ubuntu 22.04
machine and **never run on GitHub** — the file's own header says exactly what
that leaves unverified.

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
                    rc b1_nuclear cluster_config fsi triton_sf
src/core src/lhapdf src/hepmc src/pythia
python/             bindings.cpp, lipolgen/{__init__,export,cli,configs}.py, tests/
tests/              doctest suites (one file per module) + reference-table comparison
validation/         dump_polligen_reference.py, reference/*.json, VMC reconciliation
data/vmc/           ANL VMC α+d / α+t overlaps and momentum distributions (provenance in README)
docs/               DEVELOPMENT_PLAN, CONVENTIONS, USAGE, HEPMC3_CONVENTION, PYTHIA_BRIDGE,
                    T2_CHAIN, code review, OPEN_ITEMS_SOLUTIONS, PHYSICS_CHANNELS,
                    theory/ (SPIN32_FINITE_GAMMA), surveys/, open_items/ (design notes,
                    measured-number tables and gates per run)
```

## Conventions (see `docs/CONVENTIONS.md`)

Head-on frame (ion +z, electron −z); populations ordered m = +J … −J;
`TENSOR_LL_SIGN = −1` (A_zz = −(2/3) b₁/F₁, Cosyn Eq. 27 / HERMES);
⁶Li effective polarization 0.811/3 per nucleon (cluster picture);
`EmcBaseline::Epps21`; γ-matched beam energies (⁶Li 40.8 / 99.5 / 137.5,
⁷Li 40.8 / 99.5 / 117.9 GeV/u); counter-based RNG keyed by (seed, run, bunch,
event) — identical output for any thread count.

## Citing

PYTHIA 8.3 (Bierlich et al., SciPost Phys. Codebases 8 (2022)), HepMC3
(Buckley et al., CPC 260 (2021) 107310), LHAPDF6, and the physics inputs listed
in `docs/CONVENTIONS.md` and `docs/open_items/physics_literature.md`.
