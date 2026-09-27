#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""BENCHMARK_PLAN.md sec. 4 row 10 -- CLAS BONuS / EG1b database query:
the tagged spectator-spectrum SHAPE at A = 2 (tier T2).  STATUS: BLOCKED.

WHAT THE ROW WANTED.  A machine-readable BONuS (Baillie et al., PRL 108
(2012) 142001; Tkachenko et al., PRC 89 (2014) 045206) spectator-momentum
spectrum -- backward protons, p_s < 100 MeV/c, theta_pq > 100 deg, in the
deuteron REST frame -- to compare against the SHAPE of LiPolGen's deuteron-
control tagged spectator spectrum (`deuteron_channel()` -> `TaggedModel`,
n(k, cos theta_k) in the ion rest frame, averaged over the ion's M).

WHAT THE DATABASE RETURNED (queried 2026-09-23, read-only;
`data/clas_db_deuteron_query_2026-09-23.json`).  The CLAS Physics Database
(https://clas.sinp.msu.ru/cgi-bin/jlab/db.cgi) holds, for a deuteron target,
918 measurements in 19 experiments.  BONuS (E-03-012, eid 135) is ONE
measurement, E135M1: the ratio F2n/F2p against W* (Baillie 2012).  There is
no BONuS spectator-momentum spectrum and no Tkachenko 2014 entry.  EG1b
(eid 95 and 146) is inclusive (g1, A1, g1/F1, Gamma_1, ...).  So the table
this row needs is NOT machine-readable in the database.  Both papers are on
disk (refs/1110.2770.pdf, refs/1402.2477.pdf) and neither prints a p_s table
in its body.  Tkachenko et al. (PRC 89 (2014) 045206) ref. [71] puts "tables
of numerical results" -- the cos theta_pq spectra for 6 W* x 5 Q2 x 4 p_s
bins (70-85, 85-100, 100-120, 120-150 MeV/c) -- in the journal's
Supplemental Material, and that page answered HTTP 401 (APS login) on
2026-09-23; arXiv:1402.2477 has no ancillary files.  Note those tables are
R_D/S, data over BONuS's OWN spectator simulation, so even they give a p_s
SHAPE only together with that simulation.  To unblock: download the
Supplemental Material of S. Tkachenko et al. (CLAS BONuS), Phys. Rev. C 89
(2014) 045206, doi:10.1103/PhysRevC.89.045206 (docs/references/REFERENCES.md
sec. 3c), from https://journals.aps.org/prc/supplemental/10.1103/PhysRevC.89.045206
(subscription), or digitize the paper's spectator-spectrum figure with a
stated digitization error.

WHAT IT DID FIND.  The only deuteron spectator-momentum distribution in the
database is Deeps (eid 90; Klimenko et al., PRC 73 (2006) 035212,
nucl-ex/0510032): F2N x P(p_s, cos theta_pq) at p_s = 0.30-0.53 GeV/c.
VENDORED in `data/clas_deeps_klimenko2006_f2P.json` and NOT wired: it is the
high-momentum tail (280-600 MeV/c) where the paper finds PWIA adequate only
at cos theta_pq < -0.3 and FSI dominant at transverse angles, so a shape
comparison there tests LiPolGen's n(k) tail AND its FSI model together -- a
different row from this one, proposed in the record.  WHAT THE FILE HOLDS
(counted by `findings()` on every run; tabulated 2026-09-26): 115 blocks
(database measurement pages, mid 1-115), of which
  * 60 (mid 1-60) are the paper's bins: Q2 = 1.8 and 2.8 GeV^2 x six W*
    (0.94, 1.25, 1.5, 1.73, 2.02, 2.4 GeV) x five p_s -- one block each;
  * 55 (mid 61-115) carry database Q2 labels 0.18 (5 blocks) and 0.25,
    0.35, 0.45, 0.55, 0.65 (10 each), all at W* = 2.4, that are NOT the
    paper's Q2 bins (kept verbatim, NOT interpreted -- the file's own
    provenance).  Their (q2, W*, p_s) labels are not even unique: 25 label
    triples occur TWICE with different data (e.g. mid 62 and 67 are both
    (0.25, 2.4, 300 MeV), 11 and 8 rows), so nothing may key on them;
  * one block, mid 78 (labelled 0.25, 2.4, 340 MeV), has NO rows;
  * 60 rows (in 59 blocks, 55 of them a block's first cos theta_pq bin)
    read value = 0 with stat = 0 and a nonzero syst: vendored as the
    database lists them, but to be read as MISSING, not measured zeros.
A harness that wires this table must use the 60 paper-bin blocks, drop the
zero/zero rows, and state what it does with the other 55.

THE FRAME CAVEAT (02_data_nucleon_deuteron.md F-9 / sec. 5.5), which rides
with any future number from this file: BONuS tags BACKWARD in the target
rest frame at p_s < 100 MeV/c; the EIC tags a boosted spectator in the far-
forward acceptance.  Same physics, different frames; acceptance-driven
systematics do not transfer.  A shape comparison must be made in the rest
frame on both sides (LiPolGen's k IS the rest-frame relative momentum), and
BONuS's own spectrum is light-cone-corrected and acceptance-folded.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations -- NOT APPLICABLE (unpolarized, no b1).
  2. isoscalar denominator -- NOT APPLICABLE as stated; the analogue is that
     a spectator spectrum at fixed W* carries F2 of the struck NEUTRON, which
     must be divided out (or held fixed) before a shape comparison.
  3. the BONuS/Deeps observable is a CROSS SECTION on the DEUTERON (per
     nucleus); LiPolGen's n(k) is normalised to 1 per deuteron
     (integral n k^2 dk dOmega = 1).  A SHAPE comparison normalises both to
     the same integral over the common p_s window, so no absolute factor is
     compared.

INTEGRITY (since 2026-09-26).  The sha256 of every byte of both vendored
files is recorded below (`QUERY_SHA256`, `DEEPS_SHA256`, taken from the files
as committed at f0a8f1e) and re-checked on EVERY read: an altered file is
refused with a RuntimeError, so `__main__` exits 2.  The files are UTF-8 (the
catalogue's quantity names carry Greek letters) and are decoded as such
whatever the locale; the printout escapes what the terminal cannot encode,
so a non-UTF-8 stdout still gets the REPORT row and exit 0.

Run:  source env.sh && python3 validation/benchmarks/t2_bonus_spectator_shape.py

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT row -- pass, recorded fail and blocked alike, since the
row carries the verdict; 2 when the harness itself is broken (a missing or
altered vendored file, a missing dependency, any exception).  A harness is
a measurement, not a CI gate: the pytests are the gate.
"""

import collections
import hashlib
import json
import os
import sys
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
QUERY = os.path.join(HERE, "data", "clas_db_deuteron_query_2026-09-23.json")
DEEPS = os.path.join(HERE, "data", "clas_deeps_klimenko2006_f2P.json")
#: sha256 of every byte of each vendored file (recorded 2026-09-26 from the
#: files as committed at f0a8f1e); re-checked on every read.
QUERY_SHA256 = "340e04a5ae8ada1eabf4d90c2c76c82b07bda48d28c7f43cdffd881c8ee3fa11"
DEEPS_SHA256 = "1d98fd11e1b11b1e070e5504f4919b98b87dbabfe29cdca347522ea4f8737256"
#: the paper's two Q2 bins (Klimenko 2006 Sec. IV averages), GeV^2
DEEPS_PAPER_Q2 = (1.8, 2.8)

NAME = "t2_bonus_spectator_shape"
BLOCKED_REASON = (
    "BLOCKED: the CLAS Physics Database (queried 2026-09-23) holds BONuS "
    "(E-03-012, eid 135) only as ONE measurement, E135M1 = F2n/F2p vs W* "
    "(Baillie et al., PRL 108 (2012) 142001) -- no spectator-momentum "
    "spectrum, no Tkachenko 2014 entry; EG1b (eid 95, 146) is inclusive "
    "g1/A1 only.  To unblock, download the Supplemental Material (tables "
    "of numerical results, their ref. [71]) of S. Tkachenko et al. (CLAS "
    "BONuS), Phys. Rev. C 89 (2014) 045206, doi:10.1103/PhysRevC.89.045206, "
    "https://journals.aps.org/prc/supplemental/10.1103/PhysRevC.89.045206 "
    "(HTTP 401 without an APS login on 2026-09-23; arXiv:1402.2477 has no "
    "ancillary files), into validation/benchmarks/data/.  Found instead and "
    "vendored unwired: Deeps F2N x P(p_s) at 0.30-0.53 GeV/c (eid 90, "
    "Klimenko et al., PRC 73 (2006) 035212).")


def read_vendored(path, sha256):
    """The bytes of a vendored file -- RuntimeError unless their sha256 is
    the recorded one, so a moved number is refused, never silently read."""
    with open(path, "rb") as f:
        raw = f.read()
    got = hashlib.sha256(raw).hexdigest()
    if got != sha256:
        raise RuntimeError("%s: vendored file sha256 %s != recorded %s"
                           % (path, got, sha256))
    return raw


def load_query(path=QUERY):
    return json.loads(read_vendored(path, QUERY_SHA256).decode("utf-8"))


def load_deeps(path=DEEPS):
    return json.loads(read_vendored(path, DEEPS_SHA256).decode("utf-8"))


def deeps_census(d):
    """What the Deeps file holds, counted (the module docstring's census)."""
    blocks = d["blocks"]
    paper = [b for b in blocks if b["q2"] in DEEPS_PAPER_Q2]
    labels = collections.Counter((b["q2"], b["w_star"], b["ps_MeV"])
                                 for b in blocks if b["q2"] not in DEEPS_PAPER_Q2)
    zero = [(b["mid"], i) for b in blocks for i, r in enumerate(b["rows"])
            if r[1] == 0.0 and r[2] == 0.0]
    return dict(
        deeps_blocks=len(blocks),
        deeps_blocks_paper_q2=len(paper),
        deeps_blocks_unverified_q2=len(blocks) - len(paper),
        deeps_duplicate_label_triples=sum(1 for n in labels.values() if n > 1),
        deeps_empty_blocks=[b["mid"] for b in blocks if not b["rows"]],
        deeps_zero_value_zero_stat_rows=len(zero),
    )


def _printable(text):
    """`text` as the current stdout can encode it (backslash escapes for
    what it cannot), so a non-UTF-8 terminal still gets the REPORT row."""
    enc = getattr(sys.stdout, "encoding", None) or "ascii"
    return text.encode(enc, "backslashreplace").decode(enc)


def findings(query_path=QUERY, deeps_path=DEEPS):
    q = load_query(query_path)
    exps = {e["eid"]: e for e in q["experiments"]}
    d = load_deeps(deeps_path)
    return dict(
        n_experiments=len(exps),
        n_measurements=sum(e["n_measurements"] for e in exps.values()),
        bonus=exps.get(135),
        eg1b=[exps[k] for k in (95, 146) if k in exps],
        deeps=exps.get(90),
        **deeps_census(d),
    )


def run(verbose=True):
    f = findings()
    report = dict(name=NAME, generator="not evaluated",
                  reference="none machine-readable (see reason)",
                  tolerance=None, status="blocked", reason=BLOCKED_REASON,
                  findings=f)
    if verbose:
        say = lambda line: print(_printable(line))    # noqa: E731
        say("CLAS DB, deuteron target: %d experiments, %d measurements" % (
            f["n_experiments"], f["n_measurements"]))
        b = f["bonus"]
        say("  BONuS eid 135: %d measurement(s), quantities %s" % (
            b["n_measurements"], b["quantities"]))
        for e in f["eg1b"]:
            say("  EG1b eid %d: %d measurements, quantities %s" % (
                e["eid"], e["n_measurements"], e["quantities"]))
        d = f["deeps"]
        say("  Deeps eid 90: %d measurements %s; %d F2xP(p_s) blocks vendored "
            "(unwired): %d at the paper's Q2 = 1.8/2.8, %d with unverified Q2 "
            "labels (%d label triples duplicated), empty block(s) mid %s, %d "
            "value = stat = 0 rows (missing, not zeros)"
            % (d["n_measurements"], d["quantities"], f["deeps_blocks"],
               f["deeps_blocks_paper_q2"], f["deeps_blocks_unverified_q2"],
               f["deeps_duplicate_label_triples"], f["deeps_empty_blocks"],
               f["deeps_zero_value_zero_stat_rows"]))
        say("REPORT | %s | %s | %s | %s | BLOCKED" % (
            NAME, report["generator"], report["reference"], "n/a"))
        say("  " + BLOCKED_REASON)
    return report


if __name__ == "__main__":
    # Exit convention (module docstring, validation/benchmarks/README.md):
    # 0 = ran and printed its REPORT row, whatever the verdict; 2 = broken.
    try:
        run()
    except Exception:  # any exception means the harness itself is broken
        traceback.print_exc()
        sys.exit(2)
    sys.exit(0)
