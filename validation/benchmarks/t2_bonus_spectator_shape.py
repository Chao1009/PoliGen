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
nucl-ex/0510032): F2N x P(p_s, cos theta_pq) at p_s = 0.30-0.53 GeV/c,
Q2 = 1.8 and 2.8 GeV^2, six W* bins, 115 tables.  VENDORED in
`data/clas_deeps_klimenko2006_f2P.json` and NOT wired: it is the high-
momentum tail (280-600 MeV/c) where the paper finds PWIA adequate only at
cos theta_pq < -0.3 and FSI dominant at transverse angles, so a shape
comparison there tests LiPolGen's n(k) tail AND its FSI model together -- a
different row from this one, proposed in the record.

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

Run:  source env.sh && python3 validation/benchmarks/t2_bonus_spectator_shape.py

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT row -- pass, recorded fail and blocked alike, since the
row carries the verdict; 2 when the harness itself is broken (a missing or
altered vendored file, a missing dependency, any exception).  A harness is
a measurement, not a CI gate: the pytests are the gate.
"""

import json
import os
import sys
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
QUERY = os.path.join(HERE, "data", "clas_db_deuteron_query_2026-09-23.json")
DEEPS = os.path.join(HERE, "data", "clas_deeps_klimenko2006_f2P.json")

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


def findings():
    with open(QUERY) as f:
        q = json.load(f)
    exps = {e["eid"]: e for e in q["experiments"]}
    with open(DEEPS) as f:
        d = json.load(f)
    return dict(
        n_experiments=len(exps),
        n_measurements=sum(e["n_measurements"] for e in exps.values()),
        bonus=exps.get(135),
        eg1b=[exps[k] for k in (95, 146) if k in exps],
        deeps=exps.get(90),
        deeps_blocks=len(d["blocks"]),
    )


def run(verbose=True):
    f = findings()
    report = dict(name=NAME, generator="not evaluated",
                  reference="none machine-readable (see reason)",
                  tolerance=None, status="blocked", reason=BLOCKED_REASON,
                  findings=f)
    if verbose:
        print("CLAS DB, deuteron target: %d experiments, %d measurements" % (
            f["n_experiments"], f["n_measurements"]))
        b = f["bonus"]
        print("  BONuS eid 135: %d measurement(s), quantities %s" % (
            b["n_measurements"], b["quantities"]))
        for e in f["eg1b"]:
            print("  EG1b eid %d: %d measurements, quantities %s" % (
                e["eid"], e["n_measurements"], e["quantities"]))
        d = f["deeps"]
        print("  Deeps eid 90: %d measurements %s; %d F2xP(p_s) tables vendored"
              % (d["n_measurements"], d["quantities"], f["deeps_blocks"]))
        print("REPORT | %s | %s | %s | %s | BLOCKED" % (
            NAME, report["generator"], report["reference"], "n/a"))
        print("  " + BLOCKED_REASON)
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
