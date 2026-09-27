#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Vendor the T1 nucleon reference tables from the `nnpdf-data` wheel.

A DEV-TIME script, run once by a person; no harness and no test runs it, and
it is the only file under `validation/benchmarks/` that may need the network
(only when no local wheel is given).  It writes

    data/nnpdf_NMC_NC_NOTFIXED_D.json          NMC F2d            (T1, 02 F-1)
    data/nnpdf_NMC_NC_NOTFIXED.json            NMC F2d/F2p        (T1, 02 F-2)
    data/nnpdf_COMPASS15_NC_NOTFIXED_MUD.json  COMPASS g1d        (T1, 02 D-2)
    data/nnpdf_HERMES_NC_7GEV_ED.json          HERMES g1d         (T1, 02 D-1)
    data/nnpdf_E143_NC_NOTFIXED_ED.json        E143 g1d           (T1, 02 D-6)
    data/nnpdf_SMC_NC_NOTFIXED_MUD.json        SMC g1d            (recorded only)

for `t1_nmc_f2d.py`, `t1_nmc_f2d_over_f2p.py` and `t1_g1d_world.py`.

WHERE THE NUMBERS COME FROM.  The NNPDF collaboration's `nnpdf-data` package
(PyPI, GPL-3.0-or-later) ships, per data set, `metadata.yaml`, `data*.yaml`
(central values), `kinematics*.yaml` and `uncertainties*.yaml`, which NNPDF
built from the HEPData records named in each metadata file (HEPData tables
are CC0).  hepdata.net, INSPIRE and arXiv are unreachable from the machine
this was written on; PyPI is not, so the HEPData numbers arrive through that
package.  The script copies them VERBATIM (central values, kinematics, every
uncertainty component with NNPDF's ADD/MULT and CORR/UNCORR labels) into one
JSON per set with a provenance block, and computes a body sha256 that the
harnesses re-check on every read.  It transforms nothing: no unit change, no
re-ordering, no combination of uncertainties.

THE WHEEL IS PINNED.  nnpdf_data 4.1.5, py3-none-any, sha256 WHEEL_SHA256
below (the digest PyPI publishes for that file, checked 2026-09-26).  A wheel
with any other digest is refused (exit 2), and every member read is also
checked against the wheel's own RECORD.

Usage (from the repository root, after `source env.sh`):

    python validation/benchmarks/vendor_nnpdf_commondata.py \
        --wheel /path/to/nnpdf_data-4.1.5-py3-none-any.whl   # offline
    python validation/benchmarks/vendor_nnpdf_commondata.py  # downloads it
    python validation/benchmarks/vendor_nnpdf_commondata.py --check --wheel W

With no `--wheel`, `pip download --no-deps --only-binary=:all:
nnpdf_data==4.1.5` fetches the wheel into a temporary directory.  `--check`
writes nothing: it regenerates every file in memory and exits 1 unless each
is byte-identical to the committed one (the reproducibility gate for the
vendored files).  Needs PyYAML (dev time only; the harnesses read JSON).

Exit status: 0 written (or `--check` identical); 1 `--check` found a
difference; 2 anything broken (wrong wheel digest, RECORD mismatch, a count
that does not match the metadata, no PyYAML, a failed download).
"""

import argparse
import base64
import hashlib
import io
import json
import os
import subprocess
import sys
import tempfile
import zipfile

HERE = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(HERE, "data")

PACKAGE = "nnpdf_data"
VERSION = "4.1.5"
WHEEL_NAME = "nnpdf_data-4.1.5-py3-none-any.whl"
WHEEL_SHA256 = "6ed209add4d162e3b9e4801c7b538a28867e479b6ffd4781e774526a3490269a"
WHEEL_URL = ("https://files.pythonhosted.org/packages/1e/86/"
             "e8ac2b9ac3abab3e121fd4630e166e0892009e76032547b967d9e47c1e6e/"
             "nnpdf_data-4.1.5-py3-none-any.whl")
PYPI_URL = "https://pypi.org/project/nnpdf-data/4.1.5/"
WHEEL_UPLOADED = "2026-07-02T10:40:37Z"
FETCH_DATE = "2026-09-26"
COMMONDATA = "nnpdf_data/commondata/"

# ------------------------------------------------------------------ the text
# What each file says about its OWN convention.  Everything marked UNVERIFIED
# could not be read from the paper (arXiv / INSPIRE / HEPData unreachable on
# 2026-09-26); the metadata.yaml does not state it either.

_LABELS = ("ADD/MULT and CORR/UNCORR are NNPDF's treatment labels in the "
           "commondata files, copied verbatim; they are NNPDF's reading of the "
           "HEPData record, not the collaboration's words.  Uncertainties are "
           "ABSOLUTE (same units as the central value), as given; a MULT "
           "component is the relative error already multiplied by the "
           "central value.")

_G1D_DSTATE = (
    "UNVERIFIED.  metadata.yaml does not say whether the published g1d is "
    "corrected for the deuteron D state.  06_critic.md sec. 5.1 (the g1d row, "
    "line 524): 'published g1d is conventionally NOT corrected for the D "
    "state; the (1 - 1.5 omega_D) factor is applied when extracting g1n. "
    "Applying it twice, or not at all, is a 7 % error on the isoscalar "
    "combination. Verify per paper.'  The paper could not be read here, so "
    "the harness reports the tree's convention (g1_nucleus(deuteron)/2, which "
    "carries 1 - 1.5 P_D) AND the no-D-state variant (g1p + g1n)/2.")

_G1D_PER_NUCLEON = (
    "PER NUCLEON: NNPDF's y_label for the set is 'g_{1,N}(x, Q^2)' "
    "(metadata.yaml), and 06_critic.md sec. 5.1 lists g1d as per nucleon for "
    "HERMES, COMPASS and E143.  The paper's own statement is UNVERIFIED "
    "offline.")

_G1D_Q2 = ("Each point is at its own <Q^2> (the kinematics 'mid' value); no "
           "evolution to a common Q^2 is stated in the source (UNVERIFIED).")

SETS = [
    dict(
        setname="NMC_NC_NOTFIXED_D",
        out="nnpdf_NMC_NC_NOTFIXED_D.json",
        data="data_EM-F2-HEPDATA.yaml",
        kinematics="kinematics_EM-F2-HEPDATA.yaml",
        uncertainties={"EM-F2-HEPDATA": "uncertainties_EM-F2-HEPDATA.yaml"},
        primary="EM-F2-HEPDATA",
        not_vendored={},
        observable="F2 of the deuteron (F2d), per nucleon",
        survey_row="docs/benchmarking/02_data_nucleon_deuteron.md sec. 3 row "
                   "F-1 (HEPData 10.17182/hepdata.32752)",
        paper="M. Arneodo et al. (NMC), 'Measurement of the proton and "
              "deuteron structure functions, F2p and F2d, and of the ratio "
              "sigma_L/sigma_T', Nucl. Phys. B 483 (1997) 3, hep-ph/9610231 "
              "(journal and arXiv id from metadata.yaml; title as the survey "
              "row F-1 abbreviates it, not re-read)",
        convention=dict(
            normalisation="PER NUCLEON (the deuteron's average nucleon).  The "
                          "central values (0.25-0.46 at x = 0.0045-0.035) are "
                          "the size of F2 of ONE nucleon, not of two, so a "
                          "per-deuteron reading is excluded by the numbers "
                          "themselves; 06_critic.md sec. 5.1 lists NMC F2d as "
                          "per nucleon.",
            isoscalar="The deuteron is isoscalar: no non-isoscalarity "
                      "correction applies and none is stated.",
            nuclear_effects="NOT corrected for nuclear effects in the "
                            "deuteron (Fermi motion, binding, shadowing) -- "
                            "the standard for F2d, UNVERIFIED offline.",
            beam_energies="metadata.yaml version_comment: 'Hepdata "
                          "implementation of F2 deuteron averaged over "
                          "different sqrt(s)'.",
            r_in_extraction="The R = sigma_L/sigma_T used to extract F2 from "
                            "the cross section is not stated in the source "
                            "(UNVERIFIED).",
            normalisation_uncertainty="Whether NMC's overall normalisation "
                                      "uncertainty is inside 'sys' is not "
                                      "stated (UNVERIFIED).",
            uncertainty_labels=_LABELS,
        ),
    ),
    dict(
        setname="NMC_NC_NOTFIXED",
        out="nnpdf_NMC_NC_NOTFIXED.json",
        data="data.yaml",
        kinematics="kinematics.yaml",
        uncertainties={"hepdata": "uncertainties_hepdata.yaml",
                       "legacy": "uncertainties_legacy_EM-F2.yaml"},
        primary="hepdata",
        not_vendored={
            "uncertainties_EM-F2_sys_D_DEFAULT.yaml":
                "the same content, once parsed, as "
                "uncertainties_legacy_EM-F2.yaml (checked by this script at "
                "vendoring time); not duplicated.",
            "uncertainties_legacy_dw_EM-F2.yaml":
                "the legacy components (equal to 1e-12 absolute) PLUS 100 "
                "components of type "
                "'DEUTERON0..99': NNPDF's own deuteron-correction THEORY "
                "uncertainty, not data (both checked by this script at "
                "vendoring time); not vendored.",
        },
        identical_check=("uncertainties_EM-F2_sys_D_DEFAULT.yaml",
                         "uncertainties_legacy_EM-F2.yaml"),
        superset_check=("uncertainties_legacy_dw_EM-F2.yaml",
                        "uncertainties_legacy_EM-F2.yaml", "DEUTERON", 100),
        observable="the ratio F2d/F2p of per-nucleon structure functions",
        survey_row="docs/benchmarking/02_data_nucleon_deuteron.md sec. 3 row "
                   "F-2 (HEPData 10.17182/hepdata.32750)",
        paper="M. Arneodo et al. (NMC), 'Accurate measurement of F2d/F2p and "
              "Rd - Rp', Nucl. Phys. B 487 (1997) 3, hep-ex/9611022 (journal "
              "and arXiv id from metadata.yaml; title as survey row F-2 gives "
              "it, not re-read)",
        convention=dict(
            normalisation="A RATIO of PER-NUCLEON structure functions: it is "
                          "~1 at small x (0.98 at x = 0.0015) where F2n ~ F2p, "
                          "which a per-deuteron F2d would put near 2.  Apart "
                          "from that, convention-free.",
            isoscalar="The numerator is the deuteron's average nucleon, "
                      "(F2p + F2n)/2 in the absence of nuclear effects.",
            nuclear_effects="NOT corrected for nuclear effects in the "
                            "deuteron (UNVERIFIED offline).",
            r_in_extraction="metadata.yaml version_comment: 'Port of old "
                            "commondata and hepdata implementation. R and "
                            "Delta R are simultaneously determined from the "
                            "data.'",
            variants="'hepdata': stat + one 'sys' (labelled MULT, CORR) -- "
                     "the HEPData tables.  'legacy': the same stat plus five "
                     "ADD/CORR components sys_corr_1..5, NNPDF's port of its "
                     "pre-HEPData commondata (ported_from: NMCPD); their "
                     "quadrature sum matches 'hepdata' sys to <= 6.3e-4 "
                     "absolute.  The five sources are not named in the file "
                     "(UNVERIFIED).",
            uncertainty_labels=_LABELS,
        ),
    ),
    dict(
        setname="COMPASS15_NC_NOTFIXED_MUD",
        out="nnpdf_COMPASS15_NC_NOTFIXED_MUD.json",
        data="data.yaml",
        kinematics="kinematics.yaml",
        uncertainties={"default": "uncertainties.yaml"},
        primary="default",
        not_vendored={},
        observable="g1 of the deuteron (g1d), per nucleon, 6LiD target",
        survey_row="docs/benchmarking/02_data_nucleon_deuteron.md sec. 2 row "
                   "D-2 (HEPData 10.17182/hepdata.78374)",
        paper="COMPASS final g1d, Phys. Lett. B 769 (2017) 34, "
              "arXiv:1612.00620 (as survey row D-2 cites it; metadata.yaml "
              "gives only INSPIRE 1501480 / HEPData ins1501480, and the link "
              "between the two was not re-resolved offline)",
        convention=dict(
            normalisation=_G1D_PER_NUCLEON,
            d_state=_G1D_DSTATE,
            q2=_G1D_Q2,
            target="6LiD (survey 02 sec. 2 'a pointed aside'): the dilution "
                   "model is the experiment's, not the tree's.",
            uncertainty_labels=_LABELS,
        ),
    ),
    dict(
        setname="HERMES_NC_7GEV_ED",
        out="nnpdf_HERMES_NC_7GEV_ED.json",
        data="data.yaml",
        kinematics="kinematics.yaml",
        uncertainties={"default": "uncertainties.yaml"},
        primary="default",
        not_vendored={},
        observable="g1 of the deuteron (g1d), per nucleon",
        survey_row="docs/benchmarking/02_data_nucleon_deuteron.md sec. 2 row "
                   "D-1 (HEPData 10.17182/hepdata.11211)",
        paper="A. Airapetian et al. (HERMES), 'Precise determination of the "
              "spin structure function g1 of the proton, deuteron and "
              "neutron', Phys. Rev. D 75 (2007) 012007, hep-ex/0609039 (as "
              "survey row D-1 cites it; metadata.yaml gives only INSPIRE "
              "726689 / HEPData ins726689, not re-resolved offline)",
        convention=dict(
            normalisation=_G1D_PER_NUCLEON,
            d_state=_G1D_DSTATE,
            q2=_G1D_Q2 + "  The file carries an 'evol' component "
                         "('evolution systematic uncertainty'), whose meaning "
                         "for values quoted at their own <Q^2> is UNVERIFIED.",
            statistics="stat is 0 in every bin: the statistical covariance "
                       "is carried by sys_0 .. sys_14, labelled by NNPDF "
                       "'artificial correlated statistical uncertainty' "
                       "(ADD, CORR) -- a decomposition of a correlated "
                       "statistical covariance into 15 vectors.  The "
                       "per-point statistical error is sqrt(sum_k sys_k^2); "
                       "adding the vectors in quadrature per point DROPS the "
                       "bin-to-bin correlation.",
            uncertainty_labels=_LABELS,
        ),
    ),
    dict(
        setname="E143_NC_NOTFIXED_ED",
        out="nnpdf_E143_NC_NOTFIXED_ED.json",
        data="data.yaml",
        kinematics="kinematics.yaml",
        uncertainties={"default": "uncertainties.yaml"},
        primary="default",
        not_vendored={},
        observable="g1 of the deuteron (g1d), per nucleon",
        survey_row="docs/benchmarking/02_data_nucleon_deuteron.md sec. 2 row "
                   "D-6 (HEPData 10.17182/hepdata.22265)",
        paper="K. Abe et al. (E143), Phys. Rev. D 58 (1998) 112003, "
              "hep-ph/9802357 (as survey row D-6 cites it).  CONFLICT, "
              "UNRESOLVED OFFLINE: metadata.yaml's arXiv field is "
              "hep-ex/9705012, not hep-ph/9802357, beside INSPIRE 467140 / "
              "HEPData ins467140; which paper table 18 belongs to could not "
              "be checked.",
        convention=dict(
            normalisation=_G1D_PER_NUCLEON,
            d_state=_G1D_DSTATE,
            q2=_G1D_Q2,
            beam="28 points at Q^2 = 1.27-9.52 GeV^2; which E143 beam "
                 "energy (29.1 / 16.2 / 9.7 GeV) they are from is not stated "
                 "(UNVERIFIED).",
            sys_beam="'systematic uncertainty due to beam Normalization' "
                     "(MULT, CORR): 4.9 % of the central value, signed "
                     "with it, in every bin.",
            uncertainty_labels=_LABELS,
        ),
    ),
    dict(
        setname="SMC_NC_NOTFIXED_MUD",
        out="nnpdf_SMC_NC_NOTFIXED_MUD.json",
        data="data.yaml",
        kinematics="kinematics.yaml",
        uncertainties={"default": "uncertainties.yaml"},
        primary="default",
        not_vendored={},
        observable="g1 of the deuteron (g1d), per nucleon",
        survey_row="NONE: the survey (02_data_nucleon_deuteron.md) has no SMC "
                   "row; vendored as an optional, RECORDED-ONLY set",
        paper="NOT IDENTIFIED offline: metadata.yaml gives only INSPIRE "
              "471981 / HEPData ins471981 (13 points, table 7); no journal "
              "or arXiv field.",
        convention=dict(
            normalisation="PER NUCLEON on NNPDF's y_label 'g_{1,N}(x, Q^2)' "
                          "(metadata.yaml) ALONE: 06_critic.md sec. 5.1 does "
                          "not list SMC, and the paper is not identified "
                          "(UNVERIFIED).",
            d_state=_G1D_DSTATE,
            q2=_G1D_Q2,
            uncertainty_labels=_LABELS,
        ),
    ),
]


# ------------------------------------------------------------------ helpers

def die(msg):
    sys.stderr.write("vendor_nnpdf_commondata: %s\n" % msg)
    sys.exit(2)


def sha256_hex(raw):
    return hashlib.sha256(raw).hexdigest()


def body_sha256(body):
    """THE body digest the harnesses re-check: sha256 of the canonical JSON
    of `body` (sorted keys, no whitespace, ASCII, no NaN)."""
    canon = json.dumps(body, sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True, allow_nan=False)
    return sha256_hex(canon.encode("ascii"))


def fetch_wheel(tmp):
    cmd = [sys.executable, "-m", "pip", "download", "--no-deps",
           "--only-binary=:all:", "%s==%s" % (PACKAGE, VERSION), "-d", tmp]
    print("fetching: %s" % " ".join(cmd))
    res = subprocess.run(cmd, check=False)
    if res.returncode != 0:
        die("pip download failed (exit %d)" % res.returncode)
    path = os.path.join(tmp, WHEEL_NAME)
    if not os.path.isfile(path):
        die("pip download did not produce %s" % WHEEL_NAME)
    return path


def read_record(zf):
    rec = zf.read("%s-%s.dist-info/RECORD" % (PACKAGE, VERSION)).decode()
    out = {}
    for line in rec.splitlines():
        parts = line.rsplit(",", 2)
        if len(parts) == 3 and parts[1].startswith("sha256="):
            out[parts[0]] = parts[1][len("sha256="):]
    return out


def member(zf, record, path):
    """A wheel member's bytes, checked against the wheel's own RECORD."""
    raw = zf.read(path)
    want = record.get(path)
    got = base64.urlsafe_b64encode(hashlib.sha256(raw).digest()).rstrip(b"=")
    if want is None or got.decode() != want:
        die("%s: sha256 does not match the wheel's RECORD" % path)
    return raw


def _render(obj, level=0):
    """Deterministic, diff-friendly JSON: dicts one key per line, a `points`
    list one point per line, any other scalar list on one line."""
    pad = "  " * level
    if isinstance(obj, dict):
        if not obj:
            return "{}"
        items = []
        for k, v in obj.items():
            key = pad + "  " + json.dumps(k, ensure_ascii=True) + ": "
            if k == "points" and isinstance(v, list):
                rows = ",\n".join(pad + "    " + json.dumps(p, ensure_ascii=True,
                                                            allow_nan=False)
                                  for p in v)
                items.append(key + "[\n" + rows + "\n" + pad + "  ]")
            else:
                items.append(key + _render(v, level + 1))
        return "{\n" + ",\n".join(items) + "\n" + pad + "}"
    if isinstance(obj, list):
        if all(not isinstance(e, (dict, list)) for e in obj):
            return json.dumps(obj, ensure_ascii=True, allow_nan=False)
        rows = ",\n".join(pad + "  " + _render(e, level + 1) for e in obj)
        return "[\n" + rows + "\n" + pad + "]"
    return json.dumps(obj, ensure_ascii=True, allow_nan=False)


def build(spec, zf, record, yaml):
    base = COMMONDATA + spec["setname"] + "/"
    files = {}

    def load(name):
        path = base + name
        raw = member(zf, record, path)
        files[path] = sha256_hex(raw)
        return yaml.safe_load(io.BytesIO(raw))

    meta = load("metadata.yaml")
    obs = meta["implemented_observables"]
    if len(obs) != 1:
        die("%s: expected one implemented observable, got %d"
            % (spec["setname"], len(obs)))
    obs = obs[0]
    central = load(spec["data"])["data_central"]
    kin = load(spec["kinematics"])["bins"]
    variants = {}
    for label, fname in spec["uncertainties"].items():
        u = load(fname)
        variants[label] = (fname, u["definitions"], u["bins"])
    if "identical_check" in spec:
        a, b = spec["identical_check"]
        if load(a) != load(b):
            die("%s: %s and %s differ -- the not_vendored note is wrong"
                % (spec["setname"], a, b))
    if "superset_check" in spec:
        big_f, small_f, prefix, n_extra = spec["superset_check"]
        big, small = load(big_f), load(small_f)
        extra = [k for k in big["definitions"] if k not in small["definitions"]]
        ok = (len(extra) == n_extra
              and all(big["definitions"][k] == v
                      for k, v in small["definitions"].items())
              and all(str(big["definitions"][k]["type"]).startswith(prefix)
                      for k in extra)
              and all(all(abs(bb[k] - sb[k]) <= 1e-12 for k in sb)
                      for bb, sb in zip(big["bins"], small["bins"])))
        if not ok:
            die("%s: %s is not %s plus %d '%s*' components -- the "
                "not_vendored note is wrong" % (spec["setname"], big_f,
                                                small_f, n_extra, prefix))
    n = obs["ndata"]
    for what, seq in [("data", central), ("kinematics", kin)] + [
            ("uncertainties " + k, v[2]) for k, v in variants.items()]:
        if len(seq) != n:
            die("%s: %s has %d entries, metadata ndata = %d"
                % (spec["setname"], what, len(seq), n))

    kvars = list(kin[0].keys())
    points = []
    for i in range(n):
        p = {}
        for v in kvars:
            b = kin[i][v]
            if b.get("min") is not None:
                p[v + "_min"] = b["min"]
            p[v] = b["mid"]
            if b.get("max") is not None:
                p[v + "_max"] = b["max"]
        p["value"] = central[i]
        p["unc"] = {label: {k: bins[i][k] for k in defs}
                    for label, (_, defs, bins) in variants.items()}
        for label, (_, defs, bins) in variants.items():
            if set(bins[i]) != set(defs):
                die("%s: bin %d of %s has components %s, definitions %s"
                    % (spec["setname"], i, label, sorted(bins[i]), sorted(defs)))
        points.append(p)

    body = dict(
        setname=spec["setname"],
        observable=spec["observable"],
        nnpdf_observable_name=obs["observable_name"],
        nnpdf_description=obs["observable"]["description"],
        npoints=n,
        kinematics_variables={v: obs["kinematics"]["variables"][v]
                              for v in kvars},
        kinematics_note="each point's kinematics are the source's 'mid' "
                        "values; '<var>_min'/'<var>_max' are carried only "
                        "where the source gives them (it has null elsewhere)",
        primary_uncertainty_variant=spec["primary"],
        uncertainty_variants={label: dict(file=fname, definitions=defs)
                              for label, (fname, defs, _) in variants.items()},
        points=points,
    )
    prov = dict(
        what=("%s, copied VERBATIM from the nnpdf-data wheel by "
              "validation/benchmarks/vendor_nnpdf_commondata.py; nothing is "
              "transformed" % spec["observable"]),
        source=dict(
            package="nnpdf-data %s (PyPI, NNPDF Collaboration)" % VERSION,
            pypi=PYPI_URL,
            wheel_url=WHEEL_URL,
            wheel_sha256=WHEEL_SHA256,
            wheel_uploaded=WHEEL_UPLOADED,
            fetched=FETCH_DATE,
            files_read_note="every wheel member this script read (path "
                            "inside the wheel -> sha256 hex), each checked "
                            "against the wheel's RECORD; the ones not "
                            "vendored are named in 'not_vendored'",
            files_read=files,
            not_vendored=spec["not_vendored"],
        ),
        record=dict(
            hepdata=meta.get("hepdata", {}).get("url"),
            hepdata_version=meta.get("hepdata", {}).get("version"),
            hepdata_tables=obs.get("tables"),
            inspire=meta.get("iNSPIRE", {}).get("url"),
            arxiv=(meta.get("arXiv") or {}).get("url"),
            journal=(meta.get("arXiv") or {}).get("journal"),
            nnpdf_version=meta.get("version"),
            nnpdf_version_comment=meta.get("version_comment"),
            nnpdf_y_label=(obs.get("plotting") or {}).get("y_label"),
        ),
        paper=spec["paper"],
        survey_row=spec["survey_row"],
        licence=dict(
            data="CC0-1.0: HEPData records are released under CC0, and these "
                 "numbers are NNPDF's copy of the HEPData record named in "
                 "'record'",
            package="GPL-3.0-or-later (the wheel's METADATA 'License' field), "
                    "the same licence as LiPolGen",
        ),
        publishing_convention=spec["convention"],
        body_sha256_rule="sha256 of json.dumps(body, sort_keys=True, "
                         "separators=(',', ':'), ensure_ascii=True, "
                         "allow_nan=False).encode('ascii')",
    )
    doc = {"_provenance": prov, "_body_sha256": body_sha256(body),
           "body": body}
    return (_render(doc) + "\n").encode("ascii")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--wheel", help="local %s (default: pip download it)"
                    % WHEEL_NAME)
    ap.add_argument("--out-dir", default=DATA_DIR)
    ap.add_argument("--check", action="store_true",
                    help="write nothing; exit 1 unless every file is "
                         "byte-identical to the one in --out-dir")
    args = ap.parse_args(argv)
    try:
        import yaml
    except ImportError:
        die("PyYAML is needed at vendoring time (pip install pyyaml)")

    with tempfile.TemporaryDirectory() as tmp:
        wheel = args.wheel or fetch_wheel(tmp)
        raw = open(wheel, "rb").read()
        got = sha256_hex(raw)
        if got != WHEEL_SHA256:
            die("%s: sha256 %s != pinned %s" % (wheel, got, WHEEL_SHA256))
        zf = zipfile.ZipFile(io.BytesIO(raw))
        record = read_record(zf)
        differ = 0
        for spec in SETS:
            text = build(spec, zf, record, yaml)
            out = os.path.join(args.out_dir, spec["out"])
            doc = json.loads(text.decode("ascii"))
            line = ("%-40s %4d points  body sha256 %s  file sha256 %s"
                    % (spec["out"], doc["body"]["npoints"], doc["_body_sha256"],
                       sha256_hex(text)))
            if args.check:
                same = os.path.isfile(out) and open(out, "rb").read() == text
                differ += not same
                print(("IDENTICAL " if same else "DIFFERS   ") + line)
            else:
                with open(out, "wb") as f:
                    f.write(text)
                print("wrote     " + line)
    return 1 if differ else 0


if __name__ == "__main__":
    try:
        sys.exit(main())
    except SystemExit:
        raise
    except Exception:  # anything unexpected means the script is broken
        import traceback
        traceback.print_exc()
        sys.exit(2)
