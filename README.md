# LiPolGen

**Doubly polarized e + ⁶Li / e + ⁷Li DIS event generator for the EIC.**
C++17 core (spin-density kernel, cluster–spectator sampler, coherent channel),
stock **PYTHIA 8.317** as the hadronizer (through `LHAup`, no fork), **HepMC3**
output with an ion-spin attribute convention, and a **pybind11** Python module
that is a drop-in, ~100× faster replacement for the numpy generator
`PolarizedLithiumSim/evgen/polligen` it was ported from.

License: **GPL-3.0-or-later** (`LICENSE`).

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

## Status (2026-09-02)

- 275 doctest cases / 16.0 M assertions and 145 pytest cases pass.
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
Rows 5–7 and 9 are now **implemented with measured numbers** (coherent T2 via a
PYTHIA Pomeron beam, `--coherent-t2`; the Ciofi–Simula triton spectral function,
`--triton-sf`; the Glauber spectator-FSI weight, `--fsi`; the tensor-sector
radiative-correction band and radiative tails, `--rc tensor-band`); still open:
the b₁(⁶Li) convolution.

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
it is not relocatable as-is — `auditwheel repair` fixes that but pulls in
GPL-3 redistribution terms for the combined work.

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
src/core src/lhapdf src/hepmc src/pythia
python/             bindings.cpp, lipolgen/{__init__,export,cli}.py, tests/
tests/              doctest suites (one file per module) + reference-table comparison
validation/         dump_polligen_reference.py, reference/*.json, VMC reconciliation
data/vmc/           ANL VMC α+d / α+t overlaps and momentum distributions (provenance in README)
docs/               DEVELOPMENT_PLAN, CONVENTIONS, USAGE, HEPMC3_CONVENTION, PYTHIA_BRIDGE,
                    T2_CHAIN, code review, OPEN_ITEMS_SOLUTIONS, surveys/, open_items/
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
