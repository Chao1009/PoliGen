# LiPolGen

Doubly polarized e + ⁶Li / e + ⁷Li DIS event generator for the EIC.
C++17 core (spin-density kernel, cluster-spectator sampler, coherent channel),
stock PYTHIA 8.317 as hadronizer via `LHAup`, HepMC3 output, pybind11 Python
module that is a drop-in for `PolarizedLithiumSim/evgen/polligen`.

See `docs/DEVELOPMENT_PLAN.md` (feasibility verdict, architecture, work
breakdown, validation matrix), `docs/CONVENTIONS.md`, and `docs/surveys/` for
the PYTHIA 8, BeAGLE and simulation-needs surveys that motivated the design.

## Build

```bash
source env.sh                     # deps in ../deps/install
cmake -S . -B build && cmake --build build -j8
ctest --test-dir build
```

License: GPL-3.0-or-later (see LICENSE).
