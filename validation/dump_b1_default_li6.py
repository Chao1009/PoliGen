#!/usr/bin/env python3
"""Freeze `default_inclusive_kernel(LI6())`'s tensor slots as a 1e-12 gate.

Design D section 8, Agent B "Step 0".  Nothing else in the tree stores what
`pipeline.cpp::default_inclusive_kernel` returns:

  * `validation/reference/*.json` come from `dump_polligen_reference.py`,
    which dumps the PYTHON generator and builds its own kernel options with
    `toy_b1(mode="toy")` -- a different b1 from the default kernel's
    `Li6B1(MillerB1)`;
  * `tests/test_reference.cpp::build_kernel` likewise constructs its own
    `InclusiveKernel::Options` and never calls `default_inclusive_kernel`.

So the existing rtol-1e-12 gates do NOT cover that function at all, and a
regression in it -- for instance while adding the `B1Model` overload of
design D section 4.2 -- would be invisible.  This script dumps
`InclusiveKernel::tables(x, q2)`'s f1 / b1 / b2 / delta on a 200-point x grid
at three Q2 into `validation/reference/b1_default_li6.json`, which
`tests/test_b1_nuclear.cpp` (T9) reads back at rtol 1e-12.

It runs through the INSTALLED pybind11 module (`import lipolgen`), i.e. the
same C++ `default_inclusive_kernel` the test calls, so the file records what
the library does TODAY.  Regenerate it only when a default is deliberately
changed, and say so in the commit message: that is the whole point of the pin.

Floats are written with json's repr() formatter, i.e. the shortest decimal
that round-trips a double, so the C++ side (`jsonmin`, strtod) reads back the
identical bits.

Run:
    source env.sh
    python3 validation/dump_b1_default_li6.py
"""

import json
import math
import pathlib

import lipolgen as lg

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "reference"
NAME = "b1_default_li6.json"

# 200 log-spaced x from 1e-3 to 0.95 -- wide enough to cross MillerB1's own
# structure (its sign change sits near x = 0.4) and to reach the large-x tail
# where the digitized tables end.  The grid is WRITTEN INTO the file and read
# back by the test, so it never has to be reproduced by a C++ linspace.
N_X = 200
X_LO, X_HI = 1.0e-3, 0.95
Q2_GRID = [1.0, 2.5, 10.0]


def x_grid():
    lo, hi = math.log(X_LO), math.log(X_HI)
    return [math.exp(lo + (hi - lo) * i / (N_X - 1)) for i in range(N_X)]


def build():
    ion = lg.li6()
    kernel = lg.default_inclusive_kernel(ion)
    xs = x_grid()
    blocks = []
    for q2 in Q2_GRID:
        f1, b1, b2, delta = [], [], [], []
        for x in xs:
            t = kernel.tables(x, q2)
            f1.append(t.f1)
            b1.append(t.b1)
            b2.append(t.b2)
            delta.append(t.delta)
        blocks.append({"q2": q2, "f1": f1, "b1": b1, "b2": b2, "delta": delta})
    return {
        "generator": "LiPolGen/validation/dump_b1_default_li6.py",
        "description": (
            "InclusiveKernel::tables(x, q2) of default_inclusive_kernel(LI6()) "
            "with the LIBRARY DEFAULT b1 backend, B1Model::Miller = "
            "Li6B1(MillerB1) through LI6_B1_RANK2_TRANSFER, and delta_func = "
            "toy_delta_gluon(scale = 1e-2).  Design D section 8 step 0; read "
            "by tests/test_b1_nuclear.cpp T9 at rtol 1e-12.  This is the pin "
            "that makes '--b1-model miller is bit-for-bit unchanged' a "
            "testable statement."
        ),
        "kernel": {
            "function": "default_inclusive_kernel(const Ion&)",
            "ion": "6Li",
            "spin": ion.spin,
            "b1_model": "miller",
            "b1_func": "Li6B1(MillerB1)::b1_func()",
            "delta_func": "toy_delta_gluon(x, q2, f1, 1e-2)",
            "f2_source": "default (ToyF2)",
            "target_mass": kernel.target_mass,
            "tensor_gamma": kernel.tensor_gamma,
        },
        "x": x_grid(),
        "q2": Q2_GRID,
        "tables": blocks,
    }


def update_manifest(path):
    """Add this file to `_manifest.json`'s `files` list, idempotently.

    `dump_polligen_reference.py` writes the same manifest, and since the
    review of 2026-09-03 it MERGES rather than rewriting it wholesale: entries
    it did not produce (this one) and the whole `generators` map survive a
    polligen re-dump, in either order.  The `generators` map says which script
    owns which file, because the manifest's top-level `generator` key names
    only the polligen dump.
    """
    if not path.is_file():
        return
    doc = json.loads(path.read_text())
    files = doc.setdefault("files", [])
    if NAME not in files:
        files.append(NAME)
    gens = doc.setdefault("generators", {})
    gens[NAME] = "LiPolGen/validation/dump_b1_default_li6.py"
    for f in files:
        gens.setdefault(f, doc.get("generator", ""))
    path.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / NAME
    with open(out, "w") as fh:
        json.dump(build(), fh, indent=1, sort_keys=True)
        fh.write("\n")
    update_manifest(OUT / "_manifest.json")
    print("wrote %s" % out)


if __name__ == "__main__":
    main()
