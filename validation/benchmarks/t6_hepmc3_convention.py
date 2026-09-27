#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""T6 (BENCHMARK_PLAN.md sec. 2, sec. 3 "T6 not run"; 05_chain.md CH-18,
CH-19, CH-20) -- the HepMC3 files LiPolGen writes, checked CLAUSE BY CLAUSE
against the written convention, docs/HEPMC3_CONVENTION.md.

WHAT IS COMPARED.  The reference is a DOCUMENT, not data: every checkable
sentence of docs/HEPMC3_CONVENTION.md is one clause below (`CLAUSES`), each
carrying the document line(s) it checks and a short verbatim ANCHOR phrase
that must still stand on that line.  The generator side is eight small
HepMC3 Asciiv3 files written by the user's own path -- `lipolgen-run`
(`python -m lipolgen.cli`, called in-process as `lipolgen.cli.main(argv)`,
the argv recorded in the row) -- one per channel the CLI supports:
inclusive 6Li and 7Li, tagged 6Li-alpha, 7Li-alpha and d-p, coherent 6Li,
all at T2 (`--hadronize`), plus inclusive 6Li at T0 (no PYTHIA) and at T2
with `--rc tensor-band`.  Each file is read back with pyhepmc -- an
INDEPENDENT reader (pybind11 over HepMC3's ReaderAscii, not LiPolGen's
writer) -- and the raw text is read for the file-level clauses.  Every
number in a clause is computed FROM THE FILE: status codes, PDG codes, the
graph, the attributes, GenRunInfo, the four-momenta, and the charges (HepMC3
carries no charge: it is looked up per PDG code, Z from a 10LZZZAAAI code
and PYTHIA 8's own ParticleData for everything else).

NO VENDORED DATA.  The reference is the document in the tree, not a table
under data/: its sha256 is recorded on the row next to the one the clauses
were written against (DOC_SHA256_AT_AUTHORING), each clause re-finds its
anchor phrase on its cited line, and the pytest fails on a changed document
until someone re-reads the clauses against it -- so a doc change is visible.
The files are generated fresh into a temporary directory on every run.

HOW THIS DIFFERS FROM tests/test_hepmc.cpp.  That C++ suite builds ONE
synthetic tagged event by hand, writes it and reads it back with HepMC3's own
C++ ReaderAscii, and pins what the WRITER does (the round trip, the mass
rule, the RC weight names).  This harness never builds an event: it runs the
generator the way a user does, on every channel, with PYTHIA hadrons on the
record, reads the files with a different reader, and asks whether the FILES
say what the DOCUMENT says.  A divergence between the document and the code
-- in either direction -- shows up here and not there.  The two overlap on
the round-trip facts (units, masses, attribute names); only this one checks
the role table per channel, the graph example of lines 56-60, the per-event
conservation of a real hadronized record, and the document's own line
anchors.

TOLERANCE, taken from the DOCUMENT and from no file: every clause must hold
on EVERY event of EVERY file.  The numeric identities use the document's
number, 1e-9 (lines 137-138, "checks exactly this identity to `1e-9`"),
read HERE as: four-momentum, the largest |difference| over the four
components divided by max(1, summed beam energy); the Pomeron's mass^2
against its own four-vector relative to max(1, E^2); the electron mass
relative.  That four-momentum reading is NOT the check the document cites:
tests/test_hepmc.cpp compares each component with doctest's Approx,
|a - b| <= 1e-9 (1 + max(|a|, |b|)), i.e. ~1e-9 GeV absolute on px, py
(which are ~0), where the reading here allows 1e-9 x the summed beam
energy on every component -- ~6e-7 / ~7e-7 GeV on the 6Li / 7Li records
(607 / 706.5 GeV at 10 x 99.5 GeV/u), ~2e-7 GeV on d-p (210 GeV) -- i.e.
200-700x looser on the transverse components.  (Noted 2026-09-27 after a
verifier's finding, the rule NOT changed after the verdict; measured then
on the same 800 events,
not asserted: the tagged and coherent records balance to <= 4.3e-12 under
doctest's per-component rule, <= 2.2e-11 GeV absolute on E and pz and
<= 4.2e-12 GeV on px, py, and every inclusive event misses by >= 0.846
under it -- the verdict is the same under either reading.)  Charges are
exact (multiples of 1/3, compared at 1e-9).  HONESTY
NOTE: sample files were LOOKED AT before and while the checks were written
(to learn pyhepmc's API and the graph layout); no clause was dropped,
loosened or added because of what a file showed.  After the first full run
two harness slips were corrected, both visible in the record: C09's anchor
line (77 -> 78, a transcription slip) and C14's X-form test, which had
counted the unbalanced inclusive T2 events as "X carries the HFS" and so
skipped their charge -- the correction made C14 FAIL where it had silently
passed.  The remnant diagnostics were added then too; they are not clauses.
One thing seen in a file is RECORDED, not asserted: the `t` attribute is
`Event::kin.t` = pT_recoil^2 (event.hpp), which differs from the Pomeron's
|P_IP^2| = |t| + |t_min| -- the document's attribute table does not define
`t`, so this is a note, not a violation.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations -- NOT APPLICABLE: no b1 (or any structure function)
     is compared; no A-scaling assumption is used.
  2. isoscalar denominator -- NOT APPLICABLE: no nuclear PDF grid is read.
  3. per nucleus / per nucleon -- NOT APPLICABLE to the verdict: the only
     numbers compared are four-momenta and charges of ONE event record, in
     GeV in the head-on frame, per collision (e + WHOLE ion in, everything
     out).  The one cross section on the record (GenCrossSection, pb) is
     only checked to be present, finite and positive -- its per-nucleus /
     per-nucleon convention is the npz meta's business, not this document's.

WHAT THIS DOES NOT VALIDATE.  The physics of anything on the record (a
correct file of a wrong event passes); the spin_weight_<k> names in their
non-empty form (the Pipeline never fills Event::spin_weights, so the CLI
files exercise the grammar only with N = 0 -- stated in the row); the
"byte-identical with --rc off" property (a comparison between two files
across a code change; tests/test_hepmc.cpp T15/T15b own it); the producer
caveat of lines 62-68 (a rule for code, not a property of a file); the real
T6 chain legs (abconv, npsim, EICrecon, Rivet), which need containers that
are not here; and the ePIC/NuHepMC acceptance of the convention -- the
document is a PROPOSAL (its line 6), and a file can conform to a proposal
nobody downstream reads.

MEASURED 2026-09-26 (100 events per file, seed 20260926; lipolgen 0.1.0,
HepMC3 3.3.1 files, pyhepmc 2.16.1, PYTHIA 8.312): FAIL, 21 of 23 clauses
hold on all 800 events.  The two that do not are C13 and C14, the
conservation section (lines 123-138), and they fail on the four INCLUSIVE
files only (400 of 400 inclusive events): the status-1 sum misses 0.8196
(6Li) / 0.8450 (7Li) of the beam energy, and the deficit equals
P_ion - p_N -- the (A-1) remnant -- to <= 7.5e-14 on every one of them;
the missing charge is Z_ion - q(struck nucleon) on every one.  That remnant
is not written BY DESIGN (docs/USAGE.md sec. 2, "the (A-1) remnant is not
written, so the whole-nucleus balance and the total charge are deliberately
open here"; docs/CONVENTIONS.md, "Struck nucleon mass"), but THIS document
states the identity with no channel restriction, so its clause is violated
as written.  The four tagged and coherent files balance on status 1 alone
to 1.0e-13 and conserve charge exactly.  The fix belongs to the document's
owner (scope the identity, or write the remnant); nothing here moves.  4.3 s
wall for the whole row.

RE-MEASURED 2026-09-27 after decision D1 (docs/open_items/run_2026-09-27/
DECISIONS.md: the document now scopes the identity -- an inclusive record
balances per nucleon, e + struck nucleon, lines 140-147; no file changed):
PASS, 23 of 23 clauses on all 800 events, worst four-momentum residual
4.2e-13 of the beam energy.  C13/C14 were rescoped to the amended text, not
loosened: tagged and coherent files are still held to the full beams.

Run:  source env.sh && python3 validation/benchmarks/t6_hepmc3_convention.py

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT row -- pass, recorded fail and blocked alike, since the
row carries the verdict; 2 when the harness itself is broken (any
exception, lipolgen not importable included).  BLOCKED, with an
`environment:` reason, is only the HepMC3 or PYTHIA 8 tier off
(HAVE_HEPMC3, HAVE_PYTHIA8) or the pyhepmc / pythia8 module missing.  A
harness is a measurement, not a CI gate: the pytests are the gate.
"""

import contextlib
import hashlib
import io
import math
import os
import re
import sys
import tempfile
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(os.path.dirname(HERE))
DOC = os.path.join(ROOT, "docs", "HEPMC3_CONVENTION.md")
# The document these clauses were written against (2026-09-26; re-read
# 2026-09-27 after decision D1 added the per-nucleon inclusive balance,
# lines 140-147, and C13/C14 were rescoped to it).  A different sha256 is
# REPORTED on the row; the clause anchors say whether it matters.
DOC_SHA256_AT_AUTHORING = (
    "56d3a91a9bdf789b3628f9f1bc52c723887466aea3021e58ae511eb306ba4cd4")

NAME = "t6_hepmc3_convention"
N_EVENTS = 100
SEED = 20260926
TOL = 1e-9                  # the document's number, lines 137-138 (read as in the docstring)
M_E_PDG = 0.51099895e-3     # line 91
PRIMARY_STATUS = {1, 2, 3, 4}

# One file per channel the CLI supports (docs/T2_CHAIN.md, lipolgen.CHANNELS).
# `ion` is the beam's 10LZZZAAAI code, `enum` the `channel` attribute the
# writer must spell (hepmc_writer.cpp channel_name, doc line 168).
RUNS = [
    dict(label="inclusive-6Li-T2", ion=1000030060, enum="Inclusive",
         kind="inclusive", t2=True, rc=False,
         argv=["--isotope", "6Li", "--channel", "inclusive", "--hadronize"]),
    dict(label="inclusive-7Li-T2", ion=1000030070, enum="Inclusive",
         kind="inclusive", t2=True, rc=False,
         argv=["--isotope", "7Li", "--channel", "inclusive",
               "--plan", "helicity-flip", "--hadronize"]),
    dict(label="tagged-6Li-alpha-T2", ion=1000030060, enum="TaggedLi6Alpha",
         kind="alpha-tag", cluster=1000010020, t2=True, rc=False,
         argv=["--channel", "tagged-6Li-alpha", "--hadronize"]),
    dict(label="tagged-7Li-alpha-T2", ion=1000030070, enum="TaggedLi7Alpha",
         kind="alpha-tag", cluster=1000010030, t2=True, rc=False,
         argv=["--channel", "tagged-7Li-alpha", "--plan", "helicity-flip",
               "--hadronize"]),
    dict(label="tagged-d-p-T2", ion=1000010020, enum="TaggedDeuteronP",
         kind="d-tag", t2=True, rc=False,
         argv=["--channel", "tagged-d-p", "--hadronize"]),
    dict(label="coherent-6Li-T2", ion=1000030060, enum="CoherentLi6",
         kind="coherent", t2=True, rc=False,
         argv=["--channel", "coherent", "--hadronize"]),
    dict(label="inclusive-6Li-T0", ion=1000030060, enum="Inclusive",
         kind="inclusive", t2=False, rc=False,
         argv=["--isotope", "6Li", "--channel", "inclusive"]),
    dict(label="inclusive-6Li-T2-rc", ion=1000030060, enum="Inclusive",
         kind="inclusive", t2=True, rc=True,
         argv=["--isotope", "6Li", "--channel", "inclusive",
               "--rc", "tensor-band", "--hadronize"]),
]

# Event attributes of doc lines 156-179: name -> HepMC3 type.
ATTRIBUTES = [
    ("spin_J", "Double"), ("spin_M", "Double"), ("struck_cluster_m", "Double"),
    ("lam_e", "Int"), ("P_e", "Double"), ("P_z", "Double"), ("P_zz", "Double"),
    ("spin_axis_theta", "Double"), ("spin_axis_phi", "Double"),
    ("spin_category", "String"), ("run", "Int"), ("bunch", "Int"),
    ("channel", "String"), ("dis_x", "Double"), ("dis_Q2", "Double"),
    ("dis_y", "Double"), ("dis_phi", "Double"), ("spectator_k", "Double"),
    ("spectator_cos_theta", "Double"), ("spectator_phi", "Double"),
    ("alpha_s", "Double"), ("pt_s", "Double"), ("t", "Double"),
    ("x_pom", "Double"),
]

# (id, doc lines, (anchor line, verbatim ASCII phrase on it), what is checked)
CLAUSES = [
    ("C01", "16-19", (16, "Only `HepMC3::WriterAscii` (Asciiv3) is implemented"),
     "Asciiv3 only: HepMC::Version / Asciiv3-START header, END trailer, "
     "every requested event readable"),
    ("C02", "28-29", (28, "Units are GeV / mm"),
     "every event in GEV and MM"),
    ("C03", "33-37", (33, "**Primary vertex.** A single `GenVertex` with both beam"),
     "exactly two status-4 particles, e (11) and the run's ion code, both "
     "incoming at ONE vertex that has no other incoming particle"),
    ("C04", "36-37", (37, "`10LZZZAAAI` nuclear-code scheme"),
     "every nuclear PDG code is a well-formed 10LZZZAAAI (L = 0, I = 0, "
     "1 <= Z <= A)"),
    ("C05", "39-55", (55, "are hung directly off the primary vertex"),
     "every non-beam particle has a production vertex in the event; every "
     "vertex has an incoming particle"),
    ("C06", "56-60", (58, "photon and the spectator all parented by a beam particle land on the"),
     "e' and gamma* are outgoing at the primary vertex (and the alpha on an "
     "alpha tag); X is produced at its own downstream vertex fed by gamma*"),
    ("C07", "70-73", (71, "status code (`Status::Beam = 4`, `Final = 1`, `Decayed = 2`,"),
     "every status is 1, 2, 3 or 4"),
    ("C08", "72-78", (72, "`Role::HadronicX`, which is always forced"),
     "exactly one pdg-92 X per event, always status 3"),
    ("C09", "77-81", (78, "and the coherent channel's Pomeron (`Role::Pomeron`, PDG 990,"),
     "exactly one status-3 gamma* with negative (spacelike) mass; on coherent "
     "files exactly one status-3 990 with mass -sqrt|P_IP^2| from its own "
     "four-vector (1e-9), and no 990 elsewhere"),
    ("C10", "87-93", (91, "kinematics) is written with `generated_mass = 0.51099895e-3` GeV"),
     "the beam electron and e' carry generated_mass 0.51099895e-3 (rel 1e-9)"),
    ("C11", "95-109", (95, "**Role"),
     "the role table: status-3 pdg in {22, 990, 92, 2212, 2112, ions}, "
     "status-4 in {11, ions}; per channel the rows it writes (alpha tag: one "
     "status-1 alpha, one status-3 cluster ion, a status-3 nucleon; coherent: "
     "one status-1 6Li recoil; inclusive / d tag: a status-3 nucleon)"),
    ("C12", "111-121", (111, "is a genuine final-state fragment, status 1"),
     "alpha tags: the struck cluster's end vertex emits >= 1 status-1 partner "
     "(p, n or ion) and exactly one status-3 struck nucleon"),
    ("C13", "123-147", (130, "sum(status == 1 particles) + X.momentum  ==  sum(beam momenta)"),
     "four-momentum from the file: sum(status 1) == sum(beams) on a "
     "hadronized record, + X when X alone carries the HFS; never + gamma* "
     "(1e-9 of the beam energy, per component); on inclusive files the ion "
     "beam is replaced by the struck nucleon (lines 140-147, D1 2026-09-27)"),
    ("C14", "123-147 (charge)", (140, "**Inclusive channels balance per nucleon**"),
     "charge from the file's PDG codes (Z of an ion code, PYTHIA ParticleData "
     "otherwise): sum(status 1) == sum(beams) wherever X is documentation; "
     "== q_e + q_N on inclusive files"),
    ("C15", "151-179", (151, "are `GenEvent`-level attributes (id 0), named exactly"),
     "all 24 named attributes on every event, each parsing as its listed "
     "HepMC3 type"),
    ("C16", "156-168", (156, "ion spin (1 or 1.5)"),
     "spin_J in {1, 1.5}; lam_e in {-1, 0, +1}; struck_cluster_m NaN on "
     "inclusive files and finite on tagged ones; channel = the run's "
     "enumerator name"),
    ("C17", "181-189", (188, "so a reader should only ever trust index 0"),
     "a GenCrossSection on every event, index 0 finite and > 0, its error "
     "finite and >= 0"),
    ("C18", "191-193", (193, "sentinel `9.0`."),
     "a particle 'pol' attribute is a double and never the sentinel 9.0"),
    ("C19", "197-199", (197, "is built once, from the first event written, and reused"),
     "ONE GenRunInfo block in the file (one tool line, one weight-name line, "
     "both before the first event); every event reports the same names"),
    ("C20", "201-202", (201, "one `GenRunInfo::ToolInfo` entry, `{name, version}`"),
     "exactly one ToolInfo, name LiPolGen, version = lipolgen.__version__"),
    ("C21", "203-210", (203, "at index 0, plus one"),
     "weight names: 'nominal' at 0, then spin_weight_1..N consecutive; the "
     "first event carries exactly that many values"),
    ("C22", "211-228", (217, "nominal | spin_weight_1..N | rc_tensor_lo rc_tensor_hi rc_tail"),
     "--rc tensor-band: the rc names APPENDED after the spin block, exactly "
     "rc_tensor_lo rc_tensor_hi rc_tail then the _k triples, 3 (1 + N) of them"),
    ("C23", "230-232", (230, "is EMPTY when the run has `--rc off`"),
     "--rc off: no rc_* name, and 1 + N weight values"),
]


# --------------------------------------------------------------- the doc

def doc_state(path=DOC):
    """sha256 of the document and, per clause, whether its anchor phrase
    still stands on the cited line (a moved or reworded line is VISIBLE)."""
    raw = open(path, "rb").read()
    lines = raw.decode("utf-8").splitlines()
    anchors = {}
    for cid, _, (ln, phrase), _ in CLAUSES:
        ok = 1 <= ln <= len(lines) and phrase in lines[ln - 1]
        where = None
        if not ok:
            hits = [i + 1 for i, s in enumerate(lines) if phrase in s]
            where = hits[0] if hits else None
        anchors[cid] = dict(line=ln, found=ok, now_at=where)
    return dict(sha256=hashlib.sha256(raw).hexdigest(), anchors=anchors,
                n_lines=len(lines))


# ------------------------------------------------------------ PDG helpers

def nuclear_code(pid):
    """(Z, A) of a well-formed 10LZZZAAAI code with L = I = 0, else None."""
    s = str(pid)
    if len(s) != 10 or not s.startswith("10") or s[2] != "0" or s[9] != "0":
        return None
    z, a = int(s[3:6]), int(s[6:9])
    if not 1 <= z <= a:
        return None
    return z, a


def is_nuclear(pid):
    return abs(pid) >= 1000000000


class ChargeTable:
    """Charge per PDG code: Z of an ion code, PYTHIA 8's ParticleData for
    the rest.  `None` for a code neither knows (pdg 92, the string)."""

    def __init__(self):
        import pythia8                                   # noqa: WPS433
        with contextlib.redirect_stdout(io.StringIO()):
            self._py = pythia8.Pythia("", False)
        self._pd = self._py.particleData
        self.version = float(self._py.settings.parm("Pythia:versionNumber"))

    def __call__(self, pid):
        if is_nuclear(pid):
            za = nuclear_code(abs(pid))
            return None if za is None else float(za[0]) * (1 if pid > 0 else -1)
        if pid == 92 or not self._pd.isParticle(pid):
            return None
        return float(self._pd.charge(pid))


# ---------------------------------------------------------- file reading

def read_file(path):
    """(raw lines, [event dicts]) -- everything a clause needs, from the FILE."""
    import pyhepmc                                         # noqa: WPS433
    raw = open(path, encoding="utf-8", errors="replace").read().splitlines()
    events = []
    with pyhepmc.open(path) as f:
        for gev in f:
            parts = []
            for q in gev.particles:
                m = q.momentum
                pv, ev_ = q.production_vertex, q.end_vertex
                parts.append(dict(
                    id=q.id, pid=q.pid, status=q.status,
                    p4=(m.e, m.px, m.py, m.pz), m=q.generated_mass,
                    prod=pv.id if pv is not None else 0,
                    end=ev_.id if ev_ is not None else 0,
                    attrs={k: str(v) for k, v in q.attributes.items()}))
            verts = {v.id: dict(ins=[q.id for q in v.particles_in],
                                outs=[q.id for q in v.particles_out])
                     for v in gev.vertices}
            cs = gev.cross_section
            ri = gev.run_info
            events.append(dict(
                number=gev.event_number,
                units=(gev.momentum_unit.name, gev.length_unit.name),
                parts=parts, verts=verts,
                attrs={k: str(v) for k, v in gev.attributes.items()
                       if k != "GenCrossSection"},
                weights=[float(w) for w in gev.weights],
                xsec=(float(cs.xsec(0)), float(cs.xsec_err(0)))
                if cs is not None else None,
                names=list(ri.weight_names) if ri is not None else None,
                tools=[(t.name, t.version) for t in ri.tools]
                if ri is not None else None))
    return raw, events


def _add(p, q):
    return tuple(a + b for a, b in zip(p, q))


def _sum(ps):
    out = (0.0, 0.0, 0.0, 0.0)
    for p in ps:
        out = _add(out, p)
    return out


def _resid(a, b, scale):
    return max(abs(x - y) for x, y in zip(a, b)) / scale


def _primary(ev):
    beams = [p for p in ev["parts"] if p["status"] == 4]
    ends = {p["end"] for p in beams}
    return ends.pop() if len(ends) == 1 else None


def _scattered_electron(ev):
    """The status-1 e- at the primary vertex that equals k - q(gamma*): the
    graph example of line 57 says both are outgoing there."""
    v0 = _primary(ev)
    beam_e = [p for p in ev["parts"] if p["status"] == 4 and p["pid"] == 11]
    gam = [p for p in ev["parts"] if p["status"] == 3 and p["pid"] == 22]
    if v0 is None or len(beam_e) != 1 or len(gam) != 1:
        return None
    want = tuple(a - b for a, b in zip(beam_e[0]["p4"], gam[0]["p4"]))
    best, dist = None, None
    for p in ev["parts"]:
        if p["status"] == 1 and p["pid"] == 11 and p["prod"] == v0:
            d = _resid(p["p4"], want, max(1.0, beam_e[0]["p4"][0]))
            if dist is None or d < dist:
                best, dist = p, d
    return best if dist is not None and dist <= TOL else None


# ------------------------------------------------------------ the checks
#
# Each check takes (spec, raw, events, ctx) and returns (n_checked, [violation
# strings]).  n_checked = 0 means the clause does not apply to this file.

def c01(spec, raw, events, ctx):
    bad = []
    body = [s for s in raw if s.strip()]
    if not body or not body[0].startswith("HepMC::Version"):
        bad.append("line 1 is not 'HepMC::Version ...'")
    if len(body) < 2 or body[1] != "HepMC::Asciiv3-START_EVENT_LISTING":
        bad.append("line 2 is not 'HepMC::Asciiv3-START_EVENT_LISTING'")
    if not body or body[-1] != "HepMC::Asciiv3-END_EVENT_LISTING":
        bad.append("last line is not 'HepMC::Asciiv3-END_EVENT_LISTING'")
    if len(events) != ctx["n_events"]:
        bad.append("%d events read back, %d written" % (len(events), ctx["n_events"]))
    return 1, bad


def c02(spec, raw, events, ctx):
    bad = ["event %d units %s" % (e["number"], e["units"]) for e in events
           if e["units"] != ("GEV", "MM")]
    return len(events), bad


def c03(spec, raw, events, ctx):
    bad = []
    for e in events:
        beams = [p for p in e["parts"] if p["status"] == 4]
        pids = sorted(p["pid"] for p in beams)
        if pids != sorted([11, spec["ion"]]):
            bad.append("event %d status-4 pdg %s, want [11, %d]"
                       % (e["number"], pids, spec["ion"]))
            continue
        v0 = _primary(e)
        if v0 is None or v0 not in e["verts"]:
            bad.append("event %d: the beams do not end at one vertex" % e["number"])
            continue
        if sorted(e["verts"][v0]["ins"]) != sorted(p["id"] for p in beams):
            bad.append("event %d: primary vertex incoming %s"
                       % (e["number"], e["verts"][v0]["ins"]))
    return len(events), bad


def c04(spec, raw, events, ctx):
    bad, n = [], 0
    for e in events:
        for p in e["parts"]:
            if is_nuclear(p["pid"]):
                n += 1
                if p["pid"] < 0 or nuclear_code(p["pid"]) is None:
                    bad.append("event %d particle %d: pdg %d is not 10LZZZAAAI"
                               % (e["number"], p["id"], p["pid"]))
    return n, bad


def c05(spec, raw, events, ctx):
    bad = []
    for e in events:
        for p in e["parts"]:
            if p["status"] == 4:
                continue
            if p["prod"] >= 0 or p["prod"] not in e["verts"]:
                bad.append("event %d particle %d (pdg %d) has no production vertex"
                           % (e["number"], p["id"], p["pid"]))
        for vid, v in e["verts"].items():
            if not v["ins"]:
                bad.append("event %d vertex %d has no incoming particle"
                           % (e["number"], vid))
    return len(events), bad


def c06(spec, raw, events, ctx):
    bad = []
    for e in events:
        v0 = _primary(e)
        gam = [p for p in e["parts"] if p["status"] == 3 and p["pid"] == 22]
        xs = [p for p in e["parts"] if p["pid"] == 92]
        ep = _scattered_electron(e)
        if ep is None:
            bad.append("event %d: no status-1 e- at the primary vertex equal "
                       "to k - q" % e["number"])
        if len(gam) != 1 or gam[0]["prod"] != v0:
            bad.append("event %d: gamma* not outgoing at the primary vertex"
                       % e["number"])
        elif len(xs) != 1 or xs[0]["prod"] == v0 or \
                gam[0]["id"] not in e["verts"].get(xs[0]["prod"], {}).get("ins", []):
            bad.append("event %d: X not at its own vertex fed by gamma*"
                       % e["number"])
        if spec["kind"] == "alpha-tag":
            al = [p for p in e["parts"] if p["pid"] == 1000020040 and p["status"] == 1]
            if len(al) != 1 or al[0]["prod"] != v0:
                bad.append("event %d: the alpha spectator is not at the primary "
                           "vertex" % e["number"])
    return len(events), bad


def c07(spec, raw, events, ctx):
    bad = []
    for e in events:
        st = {p["status"] for p in e["parts"]}
        if not st <= PRIMARY_STATUS:
            bad.append("event %d statuses %s" % (e["number"], sorted(st)))
    return len(events), bad


def c08(spec, raw, events, ctx):
    bad = []
    for e in events:
        xs = [p for p in e["parts"] if p["pid"] == 92]
        if len(xs) != 1 or xs[0]["status"] != 3:
            bad.append("event %d: pdg-92 (status) %s"
                       % (e["number"], [p["status"] for p in xs]))
    return len(events), bad


def c09(spec, raw, events, ctx):
    bad = []
    for e in events:
        gam = [p for p in e["parts"] if p["status"] == 3 and p["pid"] == 22]
        if len(gam) != 1 or not gam[0]["m"] < 0.0:
            bad.append("event %d: status-3 gamma* %s"
                       % (e["number"], [p["m"] for p in gam]))
        pom = [p for p in e["parts"] if p["pid"] == 990]
        if spec["kind"] != "coherent":
            if pom:
                bad.append("event %d: a 990 on a non-coherent file" % e["number"])
            continue
        if len(pom) != 1 or pom[0]["status"] != 3:
            bad.append("event %d: 990 (status) %s"
                       % (e["number"], [p["status"] for p in pom]))
            continue
        E, px, py, pz = pom[0]["p4"]
        m2 = E * E - px * px - py * py - pz * pz
        mg = pom[0]["m"]
        if not (m2 < 0.0 and mg < 0.0
                and abs(mg * mg - abs(m2)) <= TOL * max(1.0, E * E)):
            bad.append("event %d: 990 mass %.12g, four-vector P^2 %.12g"
                       % (e["number"], mg, m2))
    return len(events), bad


def c10(spec, raw, events, ctx):
    bad = []
    for e in events:
        beam_e = [p for p in e["parts"] if p["status"] == 4 and p["pid"] == 11]
        ep = _scattered_electron(e)
        for tag, p in (("beam e", beam_e[0] if beam_e else None), ("e'", ep)):
            if p is None:
                bad.append("event %d: no %s" % (e["number"], tag))
            elif abs(p["m"] - M_E_PDG) > TOL * M_E_PDG:
                bad.append("event %d: %s generated_mass %.12g"
                           % (e["number"], tag, p["m"]))
    return len(events), bad


def _is_nucleon(pid):
    return pid in (2212, 2112)


def c11(spec, raw, events, ctx):
    bad = []
    for e in events:
        n = e["number"]
        for p in e["parts"]:
            if p["status"] == 3 and not (p["pid"] in (22, 990, 92, 2212, 2112)
                                         or is_nuclear(p["pid"])):
                bad.append("event %d: status-3 pdg %d not in the table" % (n, p["pid"]))
            if p["status"] == 4 and not (p["pid"] == 11 or is_nuclear(p["pid"])):
                bad.append("event %d: status-4 pdg %d not in the table" % (n, p["pid"]))
        s1 = [p["pid"] for p in e["parts"] if p["status"] == 1]
        s3 = [p["pid"] for p in e["parts"] if p["status"] == 3]
        if spec["kind"] == "alpha-tag":
            if s1.count(1000020040) != 1:
                bad.append("event %d: %d status-1 alphas" % (n, s1.count(1000020040)))
            if s3.count(spec["cluster"]) != 1:
                bad.append("event %d: %d status-3 %d clusters"
                           % (n, s3.count(spec["cluster"]), spec["cluster"]))
        if spec["kind"] == "coherent":
            if s1.count(spec["ion"]) != 1:
                bad.append("event %d: %d status-1 %d recoils" % (n, s1.count(spec["ion"]), spec["ion"]))
        elif not any(_is_nucleon(q) for q in s3):
            bad.append("event %d: no status-3 struck nucleon" % n)
    return len(events), bad


def c12(spec, raw, events, ctx):
    if spec["kind"] != "alpha-tag":
        return 0, []
    bad = []
    for e in events:
        byid = {p["id"]: p for p in e["parts"]}
        cl = [p for p in e["parts"] if p["pid"] == spec["cluster"] and p["status"] == 3]
        if len(cl) != 1 or cl[0]["end"] not in e["verts"]:
            bad.append("event %d: struck cluster has no end vertex" % e["number"])
            continue
        outs = [byid[i] for i in e["verts"][cl[0]["end"]]["outs"]]
        partners = [p for p in outs if p["status"] == 1
                    and (_is_nucleon(p["pid"]) or is_nuclear(p["pid"]))]
        struck = [p for p in outs if p["status"] == 3 and _is_nucleon(p["pid"])]
        if not partners or len(struck) != 1:
            bad.append("event %d: cluster vertex emits %s"
                       % (e["number"], [(p["pid"], p["status"]) for p in outs]))
    return len(events), bad


def _balance(e, spec=None):
    """(residual with status 1 alone, residual with status 1 + X).  On an
    inclusive file the ion beam is replaced by the struck nucleon (doc
    lines 140-147, decision D1 of 2026-09-27)."""
    beams = _initial_state(e, spec)
    fin = [p["p4"] for p in e["parts"] if p["status"] == 1]
    xs = [p["p4"] for p in e["parts"] if p["pid"] == 92]
    tin = _sum(beams)
    scale = max(1.0, tin[0])
    r1 = _resid(_sum(fin), tin, scale)
    rx = _resid(_sum(fin + xs[:1]), tin, scale) if len(xs) == 1 else float("inf")
    return r1, rx


def _initial_state(e, spec):
    """The four-vectors the record must balance against: both beams, or on
    an inclusive file the beam electron and the struck nucleon."""
    beams = [p for p in e["parts"] if p["status"] == 4]
    if spec is not None and spec["kind"] == "inclusive":
        pn = _struck_nucleon(e)
        if pn is not None:
            return [p["p4"] for p in beams if p["pid"] == 11] + [pn["p4"]]
    return [p["p4"] for p in beams]


def _x_carries_hfs(spec, e):
    """True when X is the ONLY representation of the hadronic system: a T0
    file, or a T2 event PYTHIA vetoed (the status-1 + X form balances and the
    status-1 form does not).  Otherwise X is documentation (lines 72-78)."""
    if not spec["t2"]:
        return True
    r1, rx = _balance(e, spec)
    return rx <= TOL < r1


def _struck_nucleon(e):
    """The status-3 nucleon feeding X's vertex (the per-nucleon target)."""
    xs = [p for p in e["parts"] if p["pid"] == 92]
    if len(xs) != 1:
        return None
    ins = e["verts"].get(xs[0]["prod"], {}).get("ins", [])
    byid = {p["id"]: p for p in e["parts"]}
    nuc = [byid[i] for i in ins if i in byid and byid[i]["status"] == 3
           and _is_nucleon(byid[i]["pid"])]
    return nuc[0] if len(nuc) == 1 else None


def _remnant_diagnostic(spec, e):
    """NOT a clause of this document.  For a record that does not balance:
    is the deficit exactly P_ion - p_N, the (A-1) remnant USAGE.md sec. 2
    says the inclusive record deliberately does not write?  Returns
    (deficit E / beam E, |deficit - (P_ion - p_N)| / beam E), or None."""
    beams = {p["pid"]: p["p4"] for p in e["parts"] if p["status"] == 4}
    pn = _struck_nucleon(e)
    if pn is None or spec["ion"] not in beams or 11 not in beams:
        return None
    fin = [p["p4"] for p in e["parts"] if p["status"] == 1]
    if _x_carries_hfs(spec, e):
        fin += [p["p4"] for p in e["parts"] if p["pid"] == 92]
    tin = _add(beams[11], beams[spec["ion"]])
    deficit = tuple(a - b for a, b in zip(tin, _sum(fin)))
    remnant = tuple(a - b for a, b in zip(beams[spec["ion"]], pn["p4"]))
    scale = max(1.0, tin[0])
    return deficit[0] / scale, _resid(deficit, remnant, scale)


def c13(spec, raw, events, ctx):
    bad = []
    st = ctx.setdefault("c13", dict(worst=0.0, n_final=0, n_x=0))
    diag = ctx.setdefault("remnant", {})
    for e in events:
        r1, rx = _balance(e, spec)
        if spec["t2"] and r1 <= TOL:
            st["n_final"] += 1
            st["worst"] = max(st["worst"], r1)
        elif rx <= TOL:
            st["n_x"] += 1
            st["worst"] = max(st["worst"], rx)
        else:
            bad.append("event %d: residual %.3g (status 1) / %.3g (status 1 + X)"
                       " of the beam energy"
                       % (e["number"], r1, rx))
            d = _remnant_diagnostic(spec, e)
            row = diag.setdefault(spec["label"], dict(
                n=0, n_matched=0, deficit_min=None, deficit_max=None,
                worst_match=0.0))
            row["n"] += 1
            if d is not None:
                frac, match = d
                row["deficit_min"] = frac if row["deficit_min"] is None else min(row["deficit_min"], frac)
                row["deficit_max"] = frac if row["deficit_max"] is None else max(row["deficit_max"], frac)
                row["worst_match"] = max(row["worst_match"], match)
                row["n_matched"] += int(match <= TOL)
    return len(events), bad


def c14(spec, raw, events, ctx):
    charge = ctx.get("charge")
    if charge is None:
        raise _Blocked(ctx["charge_reason"])
    bad, n = [], 0
    st = ctx.setdefault("c14", dict(n_checked=0, n_x_form=0, n_unknown=0))
    for e in events:
        if _x_carries_hfs(spec, e):
            st["n_x_form"] += 1           # X carries the HFS: no charge on file
            continue
        qs_in = [charge(p["pid"]) for p in e["parts"] if p["status"] == 4]
        pn_in = _struck_nucleon(e) if spec["kind"] == "inclusive" else None
        if pn_in is not None:            # per-nucleon balance (doc lines 140-147)
            qs_in = [charge(11), charge(pn_in["pid"])]
        qs_out = [charge(p["pid"]) for p in e["parts"] if p["status"] == 1]
        if None in qs_in or None in qs_out:
            st["n_unknown"] += 1
            bad.append("event %d: a status-1/4 pdg with no known charge" % e["number"])
            continue
        n += 1
        st["n_checked"] += 1
        if abs(sum(qs_in) - sum(qs_out)) > TOL:
            bad.append("event %d: charge in %+g, out %+g"
                       % (e["number"], sum(qs_in), sum(qs_out)))
            # diagnostic, not a clause: is the missing charge the (A-1)
            # remnant's, Z_ion - q(struck nucleon)?
            pn = _struck_nucleon(e)
            za = nuclear_code(spec["ion"])
            row = ctx.setdefault("remnant_charge", {}).setdefault(
                spec["label"], dict(n=0, n_matched=0))
            row["n"] += 1
            if pn is not None and za is not None:
                want = za[0] - charge(pn["pid"])
                row["n_matched"] += int(abs(sum(qs_in) - sum(qs_out) - want) <= TOL)
    return n, bad


_INT = re.compile(r"^[+-]?\d+$")


def _parses(value, typ):
    if typ == "Int":
        return bool(_INT.match(value))
    if typ == "Double":
        try:
            float(value)
            return True
        except ValueError:
            return False
    return isinstance(value, str)


def c15(spec, raw, events, ctx):
    bad = []
    for e in events:
        for name, typ in ATTRIBUTES:
            if name not in e["attrs"]:
                bad.append("event %d: attribute %s missing" % (e["number"], name))
            elif not _parses(e["attrs"][name], typ):
                bad.append("event %d: %s = %r is not a %sAttribute"
                           % (e["number"], name, e["attrs"][name], typ))
    return len(events), bad


def c16(spec, raw, events, ctx):
    bad = []
    for e in events:
        a, n = e["attrs"], e["number"]
        try:
            j, lam = float(a["spin_J"]), int(a["lam_e"])
            m_str = float(a["struck_cluster_m"])
        except (KeyError, ValueError) as exc:
            bad.append("event %d: %s" % (n, exc))
            continue
        if j not in (1.0, 1.5):
            bad.append("event %d: spin_J = %g" % (n, j))
        if lam not in (-1, 0, 1):
            bad.append("event %d: lam_e = %d" % (n, lam))
        if spec["kind"] == "inclusive" and not math.isnan(m_str):
            bad.append("event %d: struck_cluster_m = %g on an inclusive file" % (n, m_str))
        if spec["kind"] in ("alpha-tag", "d-tag") and not math.isfinite(m_str):
            bad.append("event %d: struck_cluster_m = %g on a tagged file" % (n, m_str))
        if a.get("channel") != spec["enum"]:
            bad.append("event %d: channel = %r, want %r" % (n, a.get("channel"), spec["enum"]))
    return len(events), bad


def c17(spec, raw, events, ctx):
    bad = []
    for e in events:
        xs = e["xsec"]
        if xs is None:
            bad.append("event %d: no GenCrossSection" % e["number"])
        elif not (math.isfinite(xs[0]) and xs[0] > 0.0
                  and math.isfinite(xs[1]) and xs[1] >= 0.0):
            bad.append("event %d: cross section %r" % (e["number"], xs))
    return len(events), bad


def c18(spec, raw, events, ctx):
    bad, n = [], 0
    for e in events:
        for p in e["parts"]:
            if "pol" in p["attrs"]:
                n += 1
                try:
                    v = float(p["attrs"]["pol"])
                except ValueError:
                    bad.append("event %d particle %d: pol %r"
                               % (e["number"], p["id"], p["attrs"]["pol"]))
                    continue
                if v == 9.0:
                    bad.append("event %d particle %d: pol is the sentinel 9"
                               % (e["number"], p["id"]))
    return n, bad


def c19(spec, raw, events, ctx):
    bad = []
    first_e = next((i for i, s in enumerate(raw) if s.startswith("E ")), len(raw))
    t_lines = [i for i, s in enumerate(raw) if s.startswith("T ")]
    w_before = [i for i, s in enumerate(raw[:first_e]) if s.startswith("W ")]
    if len(t_lines) != 1 or t_lines[0] > first_e:
        bad.append("%d tool lines (%s)" % (len(t_lines), t_lines[:3]))
    if len(w_before) != 1:
        bad.append("%d run-level weight-name lines before the first event"
                   % len(w_before))
    names0 = events[0]["names"] if events else None
    for e in events:
        if e["names"] != names0:
            bad.append("event %d: weight names %s != the first event's"
                       % (e["number"], e["names"]))
            break
    return 1, bad


def c20(spec, raw, events, ctx):
    bad = []
    tools = events[0]["tools"] if events else None
    if tools != [("LiPolGen", ctx["version"])]:
        bad.append("tool info %r, want [('LiPolGen', %r)]" % (tools, ctx["version"]))
    return 1, bad


def _spin_block(names):
    k = 0
    while 1 + k < len(names) and names[1 + k] == "spin_weight_%d" % (k + 1):
        k += 1
    return k


def c21(spec, raw, events, ctx):
    bad = []
    names = events[0]["names"] if events else None
    if not names or names[0] != "nominal":
        return 1, ["weight names %r do not start with 'nominal'" % (names,)]
    n_spin = _spin_block(names)
    ctx.setdefault("n_spin", set()).add(n_spin)
    if len(events[0]["weights"]) != len(names):
        bad.append("first event carries %d weights, %d names"
                   % (len(events[0]["weights"]), len(names)))
    for e in events:
        if len(e["weights"]) < 1:
            bad.append("event %d carries no weight" % e["number"])
    return 1, bad


def rc_names(n_spin):
    out = ["rc_tensor_lo", "rc_tensor_hi", "rc_tail"]
    for k in range(1, n_spin + 1):
        out += ["rc_tensor_lo_%d" % k, "rc_tensor_hi_%d" % k, "rc_tail_%d" % k]
    return out


def c22(spec, raw, events, ctx):
    if not spec["rc"]:
        return 0, []
    names = events[0]["names"] if events else []
    n_spin = _spin_block(names) if names and names[0] == "nominal" else 0
    tail = names[1 + n_spin:]
    bad = []
    if tail != rc_names(n_spin):
        bad.append("names after the spin block %r, want %r" % (tail, rc_names(n_spin)))
    for e in events:
        if len(e["weights"]) != 1 + n_spin + 3 * (1 + n_spin):
            bad.append("event %d carries %d weights" % (e["number"], len(e["weights"])))
            break
    return 1, bad


def c23(spec, raw, events, ctx):
    if spec["rc"]:
        return 0, []
    names = events[0]["names"] if events else []
    n_spin = _spin_block(names) if names and names[0] == "nominal" else 0
    bad = []
    if any(s.startswith("rc_") for s in names):
        bad.append("rc_* names %r on an --rc off file" % names)
    for e in events:
        if len(e["weights"]) != 1 + n_spin:
            bad.append("event %d carries %d weights, want %d"
                       % (e["number"], len(e["weights"]), 1 + n_spin))
            break
    return 1, bad


CHECKS = dict(C01=c01, C02=c02, C03=c03, C04=c04, C05=c05, C06=c06, C07=c07,
              C08=c08, C09=c09, C10=c10, C11=c11, C12=c12, C13=c13, C14=c14,
              C15=c15, C16=c16, C17=c17, C18=c18, C19=c19, C20=c20, C21=c21,
              C22=c22, C23=c23)


class _Blocked(Exception):
    """A clause that cannot be evaluated in this environment."""


# ------------------------------------------------------------- generation

def environment():
    """Why the row, or parts of it, cannot run here (None = it can).

    lipolgen itself is NOT an environment reason: a missing or broken
    lipolgen (numpy absent, an undefined symbol) is the harness broken, so
    its ImportError propagates and `__main__` exits 2, like every other
    harness (validation/benchmarks/README.md convention 6; until
    2026-09-27 it was a BLOCKED row with exit 0).  BLOCKED is kept for a
    build tier that is off and for the pyhepmc / pythia8 modules."""
    import lipolgen as lg                                  # noqa: WPS433
    out = dict(all=None, t2=None, charge=None, lg=lg)
    if not getattr(lg, "HAVE_HEPMC3", False):
        out["all"] = "environment: lipolgen built without HepMC3 (HAVE_HEPMC3 false)"
    try:
        import pyhepmc                                     # noqa: F401,WPS433
    except Exception:
        out["all"] = out["all"] or "environment: pyhepmc not importable"
    if not getattr(lg, "HAVE_PYTHIA8", False):
        out["t2"] = "environment: lipolgen built without PYTHIA 8 (HAVE_PYTHIA8 false)"
    try:
        import pythia8                                     # noqa: F401,WPS433
    except Exception:
        out["charge"] = ("environment: pythia8 Python module (the PDG charge "
                         "table) not importable")
    return out


def argv_for(spec, path, n_events=N_EVENTS, seed=SEED):
    return (list(spec["argv"]) + ["--config", "1", "--events", str(n_events),
                                 "--seed", str(seed), "--hepmc", path,
                                 "--quiet"])


def generate(spec, path, n_events=N_EVENTS, seed=SEED):
    """One file through the CLI, in process; its stdout is swallowed."""
    from lipolgen import cli                               # noqa: WPS433
    argv = argv_for(spec, path, n_events, seed)
    try:
        with contextlib.redirect_stdout(io.StringIO()):
            rc = cli.main(argv)
    except SystemExit as exc:     # a CLI refusal is a broken harness, exit 2
        raise RuntimeError("lipolgen-run %s refused: %s"
                           % (" ".join(argv), exc)) from exc
    if rc != 0:
        raise RuntimeError("lipolgen-run %s returned %r" % (" ".join(argv), rc))
    return argv


# ------------------------------------------------------------------ run

def measure(n_events=N_EVENTS, seed=SEED, out_dir=None, runs=None):
    env = environment()
    doc = doc_state()
    runs = RUNS if runs is None else runs
    res = dict(doc=doc, n_events=n_events, seed=seed, runs=[], subrows=[],
               env={k: v for k, v in env.items() if k != "lg"})
    if env["all"]:
        res["blocked"] = env["all"]
        return res
    lg = env["lg"]
    ctx = dict(n_events=n_events, version=lg.__version__)
    if env["charge"]:
        ctx["charge"], ctx["charge_reason"] = None, env["charge"]
    else:
        ctx["charge"] = ChargeTable()
        res["pythia_version"] = ctx["charge"].version
    import pyhepmc                                         # noqa: WPS433
    res["pyhepmc_version"] = pyhepmc.__version__
    res["lipolgen_version"] = lg.__version__

    tmp = None
    if out_dir is None:
        tmp = tempfile.TemporaryDirectory(prefix="t6_hepmc3_")
        out_dir = tmp.name
    per = {cid: dict(n=0, files=0, bad=[], blocked=None, failing_files=[])
           for cid, *_ in CLAUSES}
    try:
        for spec in runs:
            rrow = dict(label=spec["label"])
            if spec["t2"] and env["t2"]:
                rrow.update(status="blocked", reason=env["t2"])
                res["runs"].append(rrow)
                continue
            path = os.path.join(out_dir, spec["label"] + ".hepmc")
            t0 = time.perf_counter()
            argv = generate(spec, path, n_events, seed)
            raw, events = read_file(path)
            rrow.update(argv="lipolgen-run " + " ".join(argv[:-2] + ["<file>", "--quiet"]),
                        n_events=len(events), bytes=os.path.getsize(path),
                        hepmc_version=raw[0].split()[-1] if raw else None)
            for cid, *_ in CLAUSES:
                try:
                    n, bad = CHECKS[cid](spec, raw, events, ctx)
                except _Blocked as exc:
                    per[cid]["blocked"] = str(exc)
                    continue
                if n:
                    per[cid]["n"] += n
                    per[cid]["files"] += 1
                per[cid]["bad"] += ["%s: %s" % (spec["label"], b) for b in bad]
                if bad:
                    per[cid]["failing_files"].append(spec["label"])
            rrow["seconds"] = time.perf_counter() - t0
            rrow["status"] = "ran"
            res["runs"].append(rrow)
    finally:
        if tmp is not None:
            tmp.cleanup()

    for cid, lines, (aline, _), what in CLAUSES:
        p = per[cid]
        if p["bad"]:
            status = "fail"
        elif p["blocked"]:
            status = "blocked"
        elif p["files"] == 0:
            status = "n/a"
        else:
            status = "pass"
        res["subrows"].append(dict(
            clause=cid, doc_lines=lines, anchor_line=aline,
            anchor_found=doc["anchors"][cid]["found"], what=what,
            status=status, n_checked=p["n"], files=p["files"],
            n_violations=len(p["bad"]), first_violation=p["bad"][0] if p["bad"] else None,
            failing_files=p["failing_files"],
            reason=p["blocked"]))
    res["c13"] = ctx.get("c13")
    res["remnant"] = ctx.get("remnant", {})
    res["remnant_charge"] = ctx.get("remnant_charge", {})
    res["c14"] = ctx.get("c14")
    res["n_spin_seen"] = sorted(ctx.get("n_spin", set()))
    return res


def verdict(res):
    """pass iff every clause holds; fail on the first violation; blocked if
    nothing failed but a run or a clause could not be evaluated here."""
    if res.get("blocked"):
        return "blocked", res["blocked"]
    subs = res["subrows"]
    fails = [s for s in subs if s["status"] == "fail"]
    if fails:
        return "fail", "%s (doc lines %s): %s" % (
            fails[0]["clause"], fails[0]["doc_lines"], fails[0]["first_violation"])
    blocked = [s for s in subs if s["status"] == "blocked"] + \
              [r for r in res["runs"] if r.get("status") == "blocked"]
    if blocked:
        return "blocked", blocked[0].get("reason")
    return "pass", None


def run(verbose=True, n_events=N_EVENTS, seed=SEED, out_dir=None, runs=None):
    t0 = time.perf_counter()
    m = measure(n_events=n_events, seed=seed, out_dir=out_dir, runs=runs)
    status, why = verdict(m)
    doc = m["doc"]
    subs = m["subrows"]
    n_pass = sum(1 for s in subs if s["status"] == "pass")
    n_eval = sum(1 for s in subs if s["status"] in ("pass", "fail"))
    n_files = sum(1 for r in m["runs"] if r.get("status") == "ran")
    n_ev = sum(r.get("n_events", 0) for r in m["runs"])
    moved = [c for c, a in doc["anchors"].items() if not a["found"]]
    failing = ["%s (%s)" % (s["clause"], ", ".join(s["failing_files"]))
               for s in subs if s["status"] == "fail"]
    report = dict(
        name=NAME,
        generator=("%d/%d clauses hold on %d files, %d events (lipolgen-run, "
                   "all CLI channels, T2 + T0 + --rc)%s"
                   % (n_pass, n_eval, n_files, n_ev,
                      "; failing: " + "; ".join(failing) if failing else ""))
        if not m.get("blocked") else "not run -- %s" % m["blocked"],
        reference=("docs/HEPMC3_CONVENTION.md sha256 %s%s; %d clauses"
                   % (doc["sha256"][:16],
                      "" if doc["sha256"] == DOC_SHA256_AT_AUTHORING
                      else " (CHANGED since the clauses were written)",
                      len(CLAUSES))),
        tolerance="every clause on every event; numeric identities at the doc's 1e-9 (L137-138), "
                  "four-momentum per component relative to max(1, summed beam E)",
        status=status,
        reason=why,
        doc_sha256=doc["sha256"],
        doc_sha256_at_authoring=DOC_SHA256_AT_AUTHORING,
        doc_changed=doc["sha256"] != DOC_SHA256_AT_AUTHORING,
        anchors_moved=moved,
        subrows=subs,
        runs=m["runs"],
        config=dict(n_events=n_events, seed=seed, tol=TOL,
                    lipolgen=m.get("lipolgen_version"),
                    pyhepmc=m.get("pyhepmc_version"),
                    pythia_charge_table=m.get("pythia_version"),
                    hepmc3_file_version=next((r.get("hepmc_version") for r in m["runs"]
                                              if r.get("hepmc_version")), None)),
        conservation=m.get("c13"), charge=m.get("c14"),
        remnant_diagnostic=m.get("remnant"),
        remnant_charge_diagnostic=m.get("remnant_charge"),
        n_spin_weights_seen=m.get("n_spin_seen"),
        note=("the `t` attribute is Event::kin.t = pT_recoil^2 (event.hpp), not the "
              "Pomeron's |P_IP^2| = |t| + |t_min| its mass carries; the doc's "
              "attribute table does not define `t` -- recorded, not asserted"),
    )
    report["runtime_s"] = time.perf_counter() - t0
    if verbose:
        for r in m["runs"]:
            if r.get("status") == "ran":
                print("file  %-22s %4d events %8d bytes  %5.2f s   %s"
                      % (r["label"], r["n_events"], r["bytes"], r["seconds"], r["argv"]))
            else:
                print("file  %-22s BLOCKED: %s" % (r["label"], r.get("reason")))
        for s in subs:
            tail = ""
            if s["status"] == "fail":
                tail = " -- %d violations in %s; first: %s" % (
                    s["n_violations"], ", ".join(s["failing_files"]),
                    s["first_violation"])
            elif s["status"] == "blocked":
                tail = " -- %s" % s["reason"]
            print("  %s  L%-17s %-7s n=%-6d %s%s%s"
                  % (s["clause"], s["doc_lines"], s["status"].upper(), s["n_checked"],
                     s["what"][:70], "" if s["anchor_found"] else
                     " [ANCHOR NOT ON LINE %d]" % s["anchor_line"], tail))
        c13, c14 = m.get("c13"), m.get("c14")
        if c13:
            print("four-momentum: worst residual %.3g of the beam energy; %d events "
                  "balance on status 1 alone, %d with X (T0 or a PYTHIA veto)"
                  % (float(c13["worst"]), c13["n_final"], c13["n_x"]))
        if c14:
            print("charge: %d events checked, %d with X carrying the "
                  "HFS (no charge on the file), %d with an unknown pdg"
                  % (c14["n_checked"], c14["n_x_form"], c14["n_unknown"]))
        for label, d in sorted((m.get("remnant") or {}).items()):
            print("diagnostic (USAGE.md sec. 2, NOT a clause of this document): "
                  "%s -- %d events do not balance; deficit %.4f .. %.4f of the "
                  "beam energy, equal to P_ion - p_N (the (A-1) remnant the "
                  "inclusive record does not write) on %d of them (worst %.2g)"
                  % (label, d["n"], float(d["deficit_min"] or 0.0),
                     float(d["deficit_max"] or 0.0), d["n_matched"],
                     float(d["worst_match"])))
        for label, d in sorted((m.get("remnant_charge") or {}).items()):
            print("diagnostic (NOT a clause): %s -- %d events miss charge; "
                  "the missing charge is Z_ion - q(struck nucleon), the "
                  "remnant's, on %d of them" % (label, d["n"], d["n_matched"]))
        print("spin_weight_<k> block length seen: %s (the Pipeline never fills "
              "Event::spin_weights)" % (m.get("n_spin_seen"),))
        print("note: " + report["note"])
        print("doc sha256 %s (at authoring %s)%s"
              % (doc["sha256"], DOC_SHA256_AT_AUTHORING,
                 "" if not moved else "; anchors moved: %s" % moved))
        print("REPORT | %s | %s | %s | %s | %s" % (
            NAME, report["generator"], report["reference"], report["tolerance"],
            status.upper()))
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
