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
against the file with the tolerances tests/test_tagged.cpp applies to the
same entries (fixed 2026-09-26; until then one rtol 1e-12 was applied to
every entry, which left li7_alpha's `p2_moment_mixture_uniform` 4 % inside
its limit, about a quarter of one ULP of its summands):

  * rtol 1e-9 (`kQuadRtol`) on `norm[*].value` and
    `p2_moment_mixture_uniform`, which are themselves finite-grid
    quadratures (the README asks for ~1e-6 on them);
  * rtol 1e-12 (`kRtol`) on every other entry, with the absolute floor
    1e-300 the C++ gives `population_integrated` (`CHECK_CLOSE_AT`);
  * the C++ criterion itself, |fresh - file| <= atol + rtol * |file|, so a
    file value of exactly 0 must be matched exactly (to 1e-300 in
    `population_integrated`).

It also FAILS -- which the C++ comparison does too, and this script did not
until 2026-09-26 -- on a NaN or an infinity on either side, on an entry of
the file LiPolGen does not produce, on an entry LiPolGen produces that the
file does not have (a key, or a longer list), and on a type mismatch.  For
each channel it prints the worst ratio |fresh - file| / allowed (<= 1
passes) and, for continuity with the numbers the run records quote, the
worst plain relative difference.  Exit status 1 on any failure, 0 otherwise.

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

A13 (2026-11-03): the fourth block, `channels.li6_alpha_vmc`, is checked
after the three by `check_vmc_block` (the dilutions at the declared 1e-8).

Run:
    source env.sh
    python3 validation/repin_tagged_from_lipolgen.py     # check only
"""

import json
import math
import pathlib
import re

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
# tests/test_tagged.cpp's tolerances, entry by entry (see the docstring).
RTOL = 1e-12        # kRtol: every entry but the two quadratures below
QUAD_RTOL = 1e-9    # kQuadRtol: TaggedModel::norm and p2_moment_mixture
POP_ATOL = 1e-300   # CHECK_CLOSE_AT(..., kRtol, 1e-300): population_integrated
_QUADRATURE_PATH = re.compile(r"/norm\[\d+\]/value|/p2_moment_mixture_uniform")


def tolerance(path):
    """(rtol, atol) tests/test_tagged.cpp applies to the model entry `path`."""
    if _QUADRATURE_PATH.fullmatch(path):
        return QUAD_RTOL, 0.0
    if path.startswith("/population_integrated["):
        return RTOL, POP_ATOL
    return RTOL, 0.0


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


def _new_acc():
    return {"ratio": 0.0, "ratio_path": "", "ratio_rtol": RTOL,
            "rel": 0.0, "rel_path": "", "entries": 0, "failures": []}


def _fail(acc, path, why, uncomparable=True):
    """Record a failure.  One that has no finite |diff| (a NaN or infinity,
    a missing or extra entry, a type mismatch) also sets the worst ratio to
    inf, so the per-channel line cannot show a finite worst on a failure."""
    acc["failures"].append((path or "/", why))
    if uncomparable and acc["ratio"] != math.inf:
        acc["ratio"], acc["ratio_path"] = math.inf, path or "/"
        acc["ratio_rtol"] = tolerance(path)[0]


def worst_rel(a, b, path="", acc=None):
    """Compare LiPolGen's tree `a` with the stored tree `b`, entry by entry.

    Returns a dict: `ratio` = the worst |a - b| / (rtol |b| + atol) with
    each entry's own (rtol, atol) from tolerance() (<= 1 passes), at
    `ratio_path` with rtol `ratio_rtol`; `rel` = the worst plain |a - b| / |b|
    (|a - b| where b == 0), at `rel_path`; `entries` = numbers compared; and
    `failures` = [(path, reason)] for every entry over its tolerance,
    every NaN or infinity on either side, every missing or extra key or
    list entry and every type mismatch.  The comparison passes iff
    `failures` is empty.
    """
    if acc is None:
        acc = _new_acc()
    if isinstance(b, dict):
        if not isinstance(a, dict):
            _fail(acc, path, "LiPolGen gives %s, the file an object"
                  % type(a).__name__)
            return acc
        for k in b:
            if k in a:
                worst_rel(a[k], b[k], path + "/" + str(k), acc)
            else:
                _fail(acc, path + "/" + str(k),
                      "in the file, missing from LiPolGen's block")
        for k in a:
            if k not in b:
                _fail(acc, path + "/" + str(k),
                      "in LiPolGen's block, missing from the file")
    elif isinstance(b, list):
        if not isinstance(a, list):
            _fail(acc, path, "LiPolGen gives %s, the file a list"
                  % type(a).__name__)
            return acc
        if len(a) != len(b):
            _fail(acc, path, "%d entries from LiPolGen, %d in the file"
                  % (len(a), len(b)))
        for i in range(min(len(a), len(b))):
            worst_rel(a[i], b[i], path + "[%d]" % i, acc)
    elif isinstance(b, (int, float)) and not isinstance(b, bool):
        if isinstance(a, bool) or not isinstance(a, (int, float)):
            _fail(acc, path, "LiPolGen gives %r, the file the number %r"
                  % (a, b))
            return acc
        acc["entries"] += 1
        fa, fb = float(a), float(b)
        if not (math.isfinite(fa) and math.isfinite(fb)):
            _fail(acc, path, "non-finite: %r from LiPolGen, %r in the file"
                  % (fa, fb))
            return acc
        rtol, atol = tolerance(path)
        d = abs(fa - fb)
        allowed = rtol * abs(fb) + atol
        if allowed > 0.0:
            ratio = d / allowed
        else:
            ratio = 0.0 if d == 0.0 else math.inf
        rel = d / abs(fb) if fb else d
        if ratio > 1.0:
            _fail(acc, path, "|diff| %.3e > allowed %.3e (rtol %g, atol %g)"
                  % (d, allowed, rtol, atol), uncomparable=False)
        if ratio > acc["ratio"]:
            acc["ratio"], acc["ratio_path"], acc["ratio_rtol"] = \
                ratio, path, rtol
        if rel > acc["rel"]:
            acc["rel"], acc["rel_path"] = rel, path
    elif a != b or type(a) is not type(b):
        # strings, booleans, null: equal or a failure.
        _fail(acc, path, "%r from LiPolGen, %r in the file" % (a, b))
    return acc


def main():
    path = OUT / NAME
    doc = json.loads(path.read_text())
    print("%s provenance: %r" % (NAME, doc.get("provenance")))
    bad = []
    for key, make in CHANNELS.items():
        stored = doc["channels"][key]["model"]
        fresh = model_dump(make(), stored)
        acc = worst_rel(fresh, stored)
        print("%s: LiPolGen vs the file, %d entries: worst |diff|/allowed "
              "%.3e (at %s, rtol %g); worst relative difference %.3e (at %s); "
              "%d failure(s)"
              % (key, acc["entries"], acc["ratio"], acc["ratio_path"],
                 acc["ratio_rtol"], acc["rel"], acc["rel_path"],
                 len(acc["failures"])))
        for where, why in acc["failures"][:10]:
            print("    FAIL %s: %s" % (where, why))
        if len(acc["failures"]) > 10:
            print("    ... and %d more" % (len(acc["failures"]) - 10))
        if acc["failures"]:
            bad.append(key)
    if bad:
        raise SystemExit("LiPolGen disagrees with %s (rtol %g on norm and "
                         "p2_moment_mixture_uniform, %g elsewhere) on: %s"
                         % (NAME, QUAD_RTOL, RTOL, ", ".join(bad)))
    print("all %d model blocks agree within tests/test_tagged.cpp's "
          "tolerances (rtol %g on norm and p2_moment_mixture_uniform, %g "
          "elsewhere); nothing written." % (len(CHANNELS), QUAD_RTOL, RTOL))


# A13 (2026-11-03): the VMC alpha-d block, `channels.li6_alpha_vmc` -- a
# REDUCED block (quantum numbers, the ANL table's P_D, the two stretched-state
# dilutions) of polligen's `li6_alpha_channel(wave="vmc")`, checked against
# LiPolGen's `li6_alpha_channel(source=VmcAV18)`.  The two dilutions at the
# DECLARED relative 1e-8 (tests/test_tagged.cpp `kVmcPortRtol`; measured
# 5.389e-9 / 3.196e-9, phase_P5_port_gate.md), P_D and the quantum numbers at
# RTOL; the label is each code's own and is not compared.  Check-only.
VMC_KEY = "li6_alpha_vmc"
VMC_PORT_RTOL = 1e-8


def check_vmc_block(doc):
    """Return the list of failures of the A13 block (empty = pass)."""
    bad = []
    blk = doc["channels"].get(VMC_KEY)
    if blk is None:
        return ["channels.%s missing from the file" % VMC_KEY]
    if blk.get("port_gate_rtol") != VMC_PORT_RTOL:
        bad.append("port_gate_rtol %r in the file, %r here"
                   % (blk.get("port_gate_rtol"), VMC_PORT_RTOL))
    ch = lg.li6_alpha_channel(source=lg.ClusterWaveSource.VmcAV18)
    model = lg.TaggedModel(ch)
    rows = [(k, float(getattr(ch, k)), RTOL)
            for k in ("j_ion", "s_struck", "s_spec", "s_channel")]
    rows.append(("P_D", float(lg.VMC_P_D_LI6), RTOL))
    rows.append(("model/vector_dilution", float(model.vector_dilution()),
                 VMC_PORT_RTOL))
    rows.append(("model/tensor_dilution", float(model.tensor_dilution()),
                 VMC_PORT_RTOL))
    for path, fresh, rtol in rows:
        node = blk
        for part in path.split("/"):
            node = node[part]
        want = float(node)
        d = abs(fresh - want)
        ok = math.isfinite(fresh) and math.isfinite(want) and \
            d <= rtol * abs(want)
        print("%s/%s: LiPolGen %.17g, file %.17g, rel %.3e (rtol %g) %s"
              % (VMC_KEY, path, fresh, want, d / abs(want) if want else d,
                 rtol, "ok" if ok else "FAIL"))
        if not ok:
            bad.append(path)
    return bad


def main_with_vmc():
    """main() (the three Hulthen model blocks), then the A13 VMC block."""
    main()
    doc = json.loads((OUT / NAME).read_text())
    bad = check_vmc_block(doc)
    if bad:
        raise SystemExit("LiPolGen disagrees with %s's %s block on: %s"
                         % (NAME, VMC_KEY, ", ".join(bad)))
    print("the %s block agrees (rtol %g on the dilutions, %g elsewhere); "
          "nothing written." % (VMC_KEY, VMC_PORT_RTOL, RTOL))


if __name__ == "__main__":
    main_with_vmc()
