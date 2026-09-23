#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""CHECK-ONLY: LiPolGen's tagged `model` blocks against `tagged.json`.

CURRENT STATE (2026-09-23)
--------------------------
`validation/reference/tagged.json` is again a polligen port gate: EVERY block
is dumped from polligen by `validation/dump_polligen_reference.py`
(provenance "polligen"), which refuses to write it from a polligen without
the i^L partial-wave phase (PolarizedLithiumSim commit 1066555, 2026-09-15).
This script no longer writes anything.  It re-computes the three channels'
`model` blocks from the installed pybind11 module (`import lipolgen`, the same
C++ the tests call) on the file's own grid and sample points and checks them
against the file at rtol 1e-12 -- the same comparison tests/test_tagged.cpp
makes, from the Python side.  Exit status 1 on any disagreement.

HISTORY (2026-09-06 to 2026-09-23)
----------------------------------
On 2026-09-06 the tagged sector's S-D interference sign was found INVERTED
(docs/benchmarking/07_cw_sign_investigation.md): the partial-wave sum needs
`phi_L = i^L psi_L`, and both `TaggedModel::build_amp2` (src/core/tagged.cpp)
and polligen's `tagged._amp2_table` summed `psi_L` without it.  LiPolGen was
fixed that day with one `(-1)^floor(L/2)`; polligen was not, so this script
RE-PINNED the `li6_alpha` and `deuteron` model blocks from the fixed C++
(provenance "LiPolGen post-fix, formerly polligen", manifest generator = this
script) and the dump script carried them through.  `li7_alpha` (one L = 1
wave, a global phase) was checked, not rewritten.  polligen took the same
phase in PolarizedLithiumSim 1066555; on 2026-09-23 the carry-through was
removed and the whole file re-dumped from polligen (the two re-pinned blocks
moved by at most 2.28e-13 relative;
docs/open_items/run_2026-09-23/phase_A_port_gate.md).  The write mode was
removed at the same time: re-pinning from LiPolGen would now overwrite a live
polligen port gate.

Run:
    source env.sh
    python3 validation/repin_tagged_from_lipolgen.py     # check only
"""

import json
import pathlib

import numpy as np

import lipolgen as lg

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "reference"
NAME = "tagged.json"

# Every channel of the file; all three are checked, none is written.
CHANNELS = {
    "li6_alpha": lambda: lg.li6_alpha_channel(),
    "li7_alpha": lambda: lg.li7_alpha_channel(),
    "deuteron": lambda: lg.deuteron_channel(),
}
RTOL = 1e-12


def nearest_cell_indices(axis, points):
    """polligen's own lookup: np.clip(np.searchsorted(axis, p) - 1, 0, n-2)."""
    return np.clip(np.searchsorted(axis, points) - 1, 0, axis.size - 2)


def model_dump(channel, old_model):
    """The `model` block, on the grid and sample points the file already has."""
    model = lg.TaggedModel(channel)
    k = np.asarray(model.k)
    c = np.asarray(model.c)
    ms_ion = list(lg.m_values(channel.j_ion))
    # The k_pts / c_pts of the existing file, so the re-pin changes VALUES and
    # never the sampling.  (polligen's defaults were [0.05 .. 0.60] x
    # [-0.9 .. 0.9]; they are read back rather than retyped.)
    k_pts = np.asarray(old_model["k_pts"], dtype=float)
    c_pts = np.asarray(old_model["c_pts"], dtype=float)
    ik = nearest_cell_indices(k, k_pts)
    ic = nearest_cell_indices(c, c_pts)

    n_entries, struck_entries, p2_entries = [], [], []
    for M in ms_ion:
        n_grid = np.asarray(model.n_of_kc_table(M))
        n_entries.append({
            "M": M,
            "n_at_k_c": [[float(n_grid[i, cc]) for cc in ic] for i in ik],
        })
        struck = np.asarray(model.struck_populations(M))   # (n_mS, nk, nc)
        ms_struck = [float(v) for v in model.m_struck_values]
        struck_entries.append({
            "M": M,
            "m_s_order": ms_struck,
            "p_at_k_c": [[[float(struck[s_i, i, cc]) for cc in ic]
                          for i in ik] for s_i in range(len(ms_struck))],
        })
        p2_entries.append({"M": M, "p2_moment": float(model.p2_moment(M))})

    out = {
        "grid": {
            "k_max": 1.2, "nk": int(k.size), "nc": int(c.size),
            "k_first": float(k[0]), "k_last": float(k[-1]),
            "c_first": float(c[0]), "c_last": float(c[-1]),
        },
        "k_pts": [float(v) for v in k_pts],
        "c_pts": [float(v) for v in c_pts],
        "n_of_kc": n_entries,
        "struck_populations": struck_entries,
        "population_integrated": [
            {"M": M, "p_m_s": [float(v) for v in model.population_integrated(M)]}
            for M in ms_ion],
        "norm": [{"M": M, "value": float(model.norm(M))} for M in ms_ion],
        "vector_dilution": float(model.vector_dilution()),
        "p2_moment": p2_entries,
        "p2_moment_mixture_uniform": float(model.p2_moment_mixture(
            [1.0 / len(ms_ion)] * len(ms_ion))),
    }
    if abs(channel.s_channel - 1.0) < 1e-9:
        out["tensor_dilution"] = float(model.tensor_dilution())
    return out


def worst_rel(a, b, path="", acc=None):
    """Worst |a-b|/|b| over two identically-shaped JSON trees."""
    if acc is None:
        acc = [0.0, ""]
    if isinstance(b, dict):
        for k in b:
            worst_rel(a[k], b[k], path + "/" + str(k), acc)
    elif isinstance(b, list):
        for i, v in enumerate(b):
            worst_rel(a[i], v, path + "[%d]" % i, acc)
    elif isinstance(b, (int, float)) and not isinstance(b, bool):
        d = abs(float(a) - float(b))
        r = d / abs(float(b)) if b else d
        if r > acc[0]:
            acc[0], acc[1] = r, path
    return acc


def main():
    path = OUT / NAME
    doc = json.loads(path.read_text())
    print("%s provenance: %r" % (NAME, doc.get("provenance")))
    bad = []
    for key, make in CHANNELS.items():
        stored = doc["channels"][key]["model"]
        fresh = model_dump(make(), stored)
        r, where = worst_rel(fresh, stored)
        print("%s: LiPolGen vs the file, worst relative difference %.3e "
              "(at %s)" % (key, r, where))
        if r > RTOL:
            bad.append(key)
    if bad:
        raise SystemExit("LiPolGen disagrees with %s at rtol %g on: %s"
                         % (NAME, RTOL, ", ".join(bad)))
    print("all %d model blocks agree at rtol %g; nothing written."
          % (len(CHANNELS), RTOL))


if __name__ == "__main__":
    main()
