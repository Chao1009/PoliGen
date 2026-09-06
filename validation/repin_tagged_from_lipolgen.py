#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Re-pin `validation/reference/tagged.json`'s MODEL blocks from LiPolGen.

WHY THIS SCRIPT EXISTS INSTEAD OF `dump_polligen_reference.py`
--------------------------------------------------------------
`tagged.json` was a port gate: numbers dumped from the Python generator
`polligen` (PolarizedLithiumSim/evgen/polligen/tagged.py) that the C++ port
had to reproduce at rtol 1e-12.  On 2026-09-06 the tagged sector's S-D
interference sign was found to be INVERTED -- against Cosyn-Weiss II
Eq. (6.12), against the repository's own quadrupole sign gate, and against
the repository's own b1 sector, which applies `phi_L = i^L psi_L` explicitly
(docs/benchmarking/07_cw_sign_investigation.md).  The fix is one phase,
`(-1)^floor(L/2)`, in `TaggedModel::build_amp2` (src/core/tagged.cpp).

`polligen.tagged._amp2_table` (tagged.py:243-248) carries the SAME missing
phase.  So the reference CANNOT be regenerated from polligen: doing that
would re-bake the inverted sign and the 1e-12 gate would keep certifying the
bug.  The reference is therefore re-pinned from the FIXED C++ library, and
the tagged model block no longer tracks polligen.

WHAT MOVES AND WHAT DOES NOT
----------------------------
Only the SPIN-1 channels' `model` blocks are rewritten -- `li6_alpha` and
`deuteron` -- and only from the installed pybind11 module
(`import lipolgen`), the same C++ the tests call, exactly as
`dump_b1_default_li6.py` does.

`li7_alpha` is a single L = 1 wave: the common i is then a global phase and
the fix moves nothing (measured: `n_of_kc` and `p2_moment` identical to the
BIT across the fix).  Its model block therefore stays polligen's, byte for
byte, and the script re-checks the C++ against it at rtol 1e-12 rather than
overwriting a port gate that is still good.

Everything else in the file -- `waves`, `base`, `beam_configs`,
`boost_spectator`, `P_D_*`, the channel scalars -- is carried through
UNTOUCHED and still comes from polligen: none of it goes through
`build_amp2`.

Every other reference JSON is untouched by this script.

NOTE ON THE LAST ULP.  LiPolGen and polligen agree on the untouched
quantities to ~1e-16 relative, not to the bit -- different summation orders
in C++ and NumPy.  That is what the gate's rtol 1e-12 has always absorbed,
and it is why the 7Li cross-check below is a 1e-12 comparison and not an
equality.

Run:
    source env.sh
    python3 validation/repin_tagged_from_lipolgen.py
    python3 validation/repin_tagged_from_lipolgen.py --check   # no writes
"""

import json
import pathlib
import sys

import numpy as np

import lipolgen as lg

HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "reference"
NAME = "tagged.json"

PROVENANCE = "LiPolGen post-fix, formerly polligen"
PROVENANCE_NOTE = (
    "The li6_alpha and deuteron `model` blocks are dumped from the FIXED "
    "LiPolGen C++ library (validation/repin_tagged_from_lipolgen.py), not "
    "from polligen.  polligen's tagged._amp2_table omits the i^L "
    "partial-wave phase and therefore carries the INVERTED S-D interference "
    "sign; regenerating those blocks from it would re-bake the bug the "
    "2026-09-06 fix removed.  See "
    "docs/benchmarking/07_cw_sign_investigation.md and "
    "docs/open_items/run_2026-09-06/phase_CW_numbers.md.  Everything else in "
    "this file -- including li7_alpha's model block, one L = 1 wave with "
    "nothing to interfere with -- is still polligen's."
)

# The spin-1 channels: S + D, so `build_amp2`'s (-1)^floor(L/2) is a REAL
# relative phase and these are the blocks the 2026-09-06 fix moves.
REPINNED = {
    "li6_alpha": lambda: lg.li6_alpha_channel(),
    "deuteron": lambda: lg.deuteron_channel(),
}
# One L = 1 wave: nothing to interfere with, nothing moves.  Checked, kept.
UNCHANGED = {
    "li7_alpha": lambda: lg.li7_alpha_channel(),
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


def update_manifest(path):
    """Point `tagged.json`'s generator at THIS script, idempotently."""
    if not path.is_file():
        return
    doc = json.loads(path.read_text())
    gens = doc.setdefault("generators", {})
    gens[NAME] = "LiPolGen/validation/repin_tagged_from_lipolgen.py"
    for f in doc.get("files", []):
        gens.setdefault(f, doc.get("generator", ""))
    path.write_text(json.dumps(doc, indent=1, sort_keys=True) + "\n")


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
    check_only = "--check" in sys.argv[1:]
    path = OUT / NAME
    doc = json.loads(path.read_text())

    # 7Li: the port gate against polligen is untouched by the fix.  Verify it
    # still holds and leave the stored block exactly where it is.
    for key, make in UNCHANGED.items():
        stored = doc["channels"][key]["model"]
        fresh = model_dump(make(), stored)
        r, where = worst_rel(fresh, stored)
        print("%s: kept polligen's block; LiPolGen agrees to %.3e (worst at %s)"
              % (key, r, where))
        if r > RTOL:
            raise SystemExit(
                "%s has ONE partial wave, so the i^L phase is global and this "
                "block must still track polligen at rtol %g.  It does not "
                "(%.3e at %s) -- stop and find out why before writing "
                "anything." % (key, RTOL, r, where))

    moved = []
    for key, make in REPINNED.items():
        old = doc["channels"][key]["model"]
        new = model_dump(make(), old)
        r, where = worst_rel(new, old)
        print("%s: re-pinned from LiPolGen; worst move vs the old file %.3e "
              "(at %s)" % (key, r, where))
        moved.append(key)
        doc["channels"][key]["model"] = new

    doc["provenance"] = PROVENANCE
    doc["provenance_note"] = PROVENANCE_NOTE
    if check_only:
        print("--check: nothing written.  Would re-pin: %s" % ", ".join(moved))
        return
    with open(path, "w") as fh:
        json.dump(doc, fh, indent=1, sort_keys=True)
        fh.write("\n")
    update_manifest(OUT / "_manifest.json")
    print("rewrote %s (model blocks re-pinned: %s)" % (path, ", ".join(moved)))


if __name__ == "__main__":
    main()
