#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""BENCHMARK_PLAN.md sec. 4 row 3 -- HERMES b1^d Table II as an assertion
(tier T2).

WHAT IS COMPARED.  The six (<x>, <Q^2>) rows of HERMES PRL 95 (2005) 242001
Table II (`data/hermes2005_b1d_table2.json`) against LiPolGen's deuteron b1,
per nucleon, at each row's own (<x>, <Q^2>), in five configurations and one
null baseline, with NO tuning of any of them:

  cdks_conv_mstw     DeuteronConvolutionB1 with unpol = MstwSF() -- the A = 2
                     gate's passing configuration (AV18 `fdeut.av18`, CDKS
                     Eq. (21) finite-|q| delta, r1998, CDKS Eq. (22) target
                     mass): the MstwSF path.  Needs the PYTHIA 8 tier.
  cdks_conv_default  DeuteronConvolutionB1() -- the class's shipped default
                     Options (unpol = ToyF2, otherwise as above).
  miller_shipped     MillerB1() -- the deuteron b1 the generator ships
                     (`--b1-model miller`): Miller's per-DEUTERON curve times
                     B1_MILLER_TABLE_TO_PER_NUCLEON = 0.5 (registry row 2,
                     status quo).
  miller_x1          the same curve WITHOUT the 0.5 (registry row 2's other
                     option), i.e. MillerB1() / B1_MILLER_TABLE_TO_PER_NUCLEON.
                     A PRICE, not a proposal.
  cdks_digitized     CdksB1() -- the digitized CDKS Fig. 4 column (Q^2 = 2.5,
                     per nucleon, B1_CDKS_TABLE_TO_PER_NUCLEON = 1), read at
                     every row's x regardless of its Q^2.
  null_b1_zero       b1 = 0 everywhere: NOT a LiPolGen model, the baseline a
                     chi2 must beat to say the data prefer a model to nothing.

A_zz^d is formed the way HERMES formed it, inverted: A_zz = -(2/3) b1/F1^d
with F1^d = (1 + gamma^2) F2^d/(2x(1 + R1990)), gamma^2 = 4 M^2 x^2/Q^2,
F2^d = (F2^p + F2^n)/2 per nucleon.  HERMES took F2^p from ALLM97 and
F2^n/F2^p from NMC 1992, neither of which is in the tree; this harness uses
MstwSF (ToyF2 if the PYTHIA tier is absent) and PRICES the difference against
the F1^d HERMES's own table implies, F1_impl = -(2/3) b1/A_zz (rounded
3-digit entries, so that ratio is itself good to ~1-5 %).  R1990 is
`data/whitlow1990_r1990.json`; the R1998 alternative is priced alongside.

TOLERANCE.  Per bin, the quadrature sum sqrt(stat^2 + syst^2) of the row
(the systematic includes a ~1e-3 A_zz normalisation that HERMES says is
correlated over bins; treated here as uncorrelated, which is what the plan
row asks).  Per configuration: chi2 = sum (model - data)^2/sigma^2 over the
six bins, ndf = 6 (no parameter is fitted).  STATUS: "pass" iff every bin is
within its tolerance; the chi2 is reported either way.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations, both sides named.  HERMES: PER NUCLEON (their F2^d =
     (F2^p + F2^n)/2).  DeuteronConvolutionB1 / CdksB1: PER NUCLEON, the 0.5
     dropped (CDKS Eq. (10)'s explicit 1/A; B1_CDKS_TABLE_TO_PER_NUCLEON = 1,
     CERTAIN).  MillerB1: Miller's table is PER DEUTERON and the shipped
     value keeps the 0.5 (B1_MILLER_TABLE_TO_PER_NUCLEON; LIKELY, NOT CERTAIN
     -- registry row 2).  A-scaling: none -- A = 2 on both sides; no 2/6
     (that is the 6Li inert-alpha assumption and is NOT applied here).  x is
     Bjorken x on the NUCLEON mass on both sides (runs to 2 for the deuteron;
     CDKS's convention, and HERMES's).
  2. isoscalar denominator: F2^d, F1^d are the isoscalar nucleon
     (F2^p + F2^n)/2 on both sides -- asserted in `measure()` by recomputing
     F1^d from f2p and f2n explicitly, never from a nuclear grid.
  3. structure functions per nucleon (no cross section enters).

Born level: HERMES's A_zz is RC-unfolded (RADGEN); the comparison is to a
Born structure function.  Q^2 = 0.51 at bin 1 is below MSTW2008's 1 GeV^2
grid floor (PYTHIA's reader extrapolates) and below R1990's SLAC x range.

INTEGRITY (since 2026-09-26).  The sha256 of every byte of each vendored
file is recorded below (`DATA_SHA256`, `R1990_SHA256`, taken from the files
as committed at f0a8f1e) and re-checked on EVERY read: an altered table is
refused with a RuntimeError, so `__main__` exits 2, never 0 on moved numbers.

Run:  source env.sh && python3 validation/benchmarks/t2_hermes_b1_table2.py [--plot]

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT row -- pass, recorded fail and blocked alike, since the
row carries the verdict; 2 when the harness itself is broken (a missing or
altered vendored file, a missing dependency, any exception).  A harness is
a measurement, not a CI gate: the pytests are the gate.
"""

import hashlib
import json
import math
import os
import sys
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "hermes2005_b1d_table2.json")
R1990_DATA = os.path.join(HERE, "data", "whitlow1990_r1990.json")
#: sha256 of every byte of each vendored file (recorded 2026-09-26 from the
#: files as committed at f0a8f1e); re-checked on every read.
DATA_SHA256 = "20b9104edf52414822e88305f7144566568ee2979bd531f77dd79f3a92006fc5"
R1990_SHA256 = "ce330eeef6305315a78d66422ff81fe8307bf71fdca4b8ef6e781e71698591d9"

NAME = "t2_hermes_b1_table2"
M_NUCLEON = 0.938272   # GeV; gamma^2 = 4 M^2 x^2 / Q^2 = Q^2/nu^2
CONFIGS = ("cdks_conv_mstw", "cdks_conv_default", "miller_shipped",
           "miller_x1", "cdks_digitized", "null_b1_zero")


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


def load_table(path=DATA):
    d = json.loads(read_vendored(path, DATA_SHA256).decode("utf-8"))
    s = d["scale"]
    out = []
    for r in d["rows"]:
        x, q2, a, sa, ya, b, sb, yb = r
        out.append(dict(x=x, q2=q2, azz=a * s, azz_err=math.hypot(sa, ya) * s,
                        b1=b * s, b1_err=math.hypot(sb, yb) * s,
                        f1_implied=-(2.0 / 3.0) * b / a))
    return out


def r1990_factory(path=R1990_DATA):
    p = json.loads(read_vendored(path, R1990_SHA256).decode("utf-8"))

    def r1990(x, q2):
        th = 1.0 + 12.0 * (q2 / (q2 + 1.0)) * (0.125 ** 2 / (0.125 ** 2 + x * x))
        return (p["b1"] / math.log(q2 / p["Lambda2_GeV2"]) * th + p["b2"] / q2
                + p["b3"] / (q2 * q2 + 0.09))
    return r1990


def have_mstw():
    import lipolgen._lipolgen as _l
    if not getattr(_l, "HAVE_PYTHIA8", False):
        return False, "PYTHIA 8 tier not built"
    d = _l.pythia8_pdfdata_dir()
    if not (d and os.path.isfile(os.path.join(d, "mstw2008lo.00.dat"))):
        return False, "PYTHIA 8 built without pdfdata/mstw2008lo.00.dat"
    return True, ""


def build_models():
    import lipolgen._lipolgen as _l
    models, blocked = {}, {}
    ok, why = have_mstw()
    if ok:
        o = _l.DeuteronConvolutionB1.Options()
        o.unpol = _l.MstwSF()
        models["cdks_conv_mstw"] = _l.DeuteronConvolutionB1(o)
    else:
        blocked["cdks_conv_mstw"] = why
    models["cdks_conv_default"] = _l.DeuteronConvolutionB1()
    mil = _l.MillerB1()
    models["miller_shipped"] = mil
    models["miller_x1"] = mil          # divided below, at evaluation
    models["cdks_digitized"] = _l.CdksB1()
    models["null_b1_zero"] = None
    f2 = _l.MstwSF() if ok else _l.ToyF2()
    return models, blocked, f2, ("MstwSF" if ok else "ToyF2")


# A5 NOTE (decision A5 of 2026-11-03 and the supervisor's note of 2026-11-06
# under A5, docs/open_items/run_2026-11-03/DECISIONS.md -- the supervising
# session's choices under the author's delegation; applied 2026-11-06, phase
# P3; validation/benchmarks/README.md convention 8).  The row's HEADLINE is
# the shipped b1, `miller_shipped` (HEADLINE; `cdks_conv_mstw` until A5, or
# `cdks_conv_default` without MSTW), and because the rule names the whole
# shipped configuration its A_zz is formed on the shipped F1^d backend and R
# as `lipolgen.make_config()` builds them: F2^d per nucleon from the shipped
# deuteron kernel (`shipped_inclusive_kernel(deuteron())`, CT18NLO at LO,
# frozen below its Q2 grid as shipped) and R = `r_sigma_lt` (the kernel's
# own; its f1 equals F2/(2x(1 + r_sigma_lt)) at every bin), in HERMES's
# construction F1 = (1 + gamma^2) F2/(2x(1 + R)).  Every other configuration,
# the declared `cdks_conv_mstw` among them, keeps its A_zz on MstwSF and
# R1990 (unchanged), recorded.  The b1 chi2 does not depend on F1^d.  Status
# FAIL -> FAIL.  RECORDED, not acted on: the shipped inclusive deuteron
# kernel's own b1 was (until J1, 2026-11-10 -- the NOTE below) 0.32753 x MillerB1() (= (2/6) x LI6_B1_RANK2_TRANSFER,
# the 6Li mapping), not MillerB1() -- `miller_shipped` is the MillerB1() the
# `--b1-model miller` backend names, as declared above.
#
# J1 NOTE (decision J1 of 2026-11-10, THE AUTHOR'S, docs/open_items/run_2026-11-10/
# DECISIONS.md; record phase_P2_deuteron_b1.md).  The inclusive deuteron kernel's
# rank-2 slot is now MillerB1() itself, per nucleon, unmapped, so the headline
# `miller_shipped` IS the b1 a default deuteron run carries and convention 8's
# default-run clause applies to this row again.  Measured 2026-11-10 on
# `shipped_inclusive_kernel(deuteron()).tables(x, q2).b1` at the six bins: equal to
# MillerB1() bin for bin, chi2(b1)/6 = 5.2619 (4/6) and A_zz 7.3080/6 (3/6) on the
# shipped F1^d -- the headline's own numbers.  The 0.32753 x MillerB1() reading
# (chi2(b1)/6 = 14.56, 2/6; A_zz 15.68/6 on the shipped F1^d) is history.  This
# harness never read the kernel's b1, so no number of this row moved.
#
# J4 NOTE (decision J4 of 2026-11-10, docs/open_items/run_2026-11-10/DECISIONS.md --
# the SUPERVISOR's choice under the author's delegation, not the author's words;
# applied 2026-11-12, record phase_P3_harness_rows.md).  The headline keeps all six
# bins, its rule and its status.  Beside it, ONE RECORDED sub-row (INGRID_SUBROW):
# the headline model's chi2(b1) and chi2(A_zz) (A_zz on the shipped F1^d, as the
# headline) on the bins whose <Q2> is at or above the shipped unpolarised grid's
# OWN Q2 floor -- read from the process's one shared CT18NLO object, the one
# `shipped_inclusive_kernel` hands the kernel (`make_lhapdf_unpol_sf(
# UNPOL_SF_DEFAULT).q2_min`), never typed.  CHOSEN AFTER THE VALUES WERE SEEN: a
# diagnostic, no verdict.  Below the floor D7 freezes F2 at the edge (the
# author's decision, not re-opened): the shipped F2^d there equals its value at
# the floor exactly (checked per bin, `f2d_frozen_at_floor`), so the b1 <-> A_zz
# conversion uses an F1^d that is not the measured one.  Measured 2026-11-12:
# floor 1.677025; below it the bins <Q2> = 0.51, 1.06, 1.65 (x = 0.012, 0.032,
# 0.063), shipped F1^d / the F1^d HERMES's own b1 and A_zz columns imply = 1.761 /
# 1.413 / 1.338; on the three in-grid bins chi2(b1)/3 = 1.894 (2/3 within) and
# chi2(A_zz)/3 = 1.178 (3/3 within).
#
# K1 NOTE (decision K1 of 2026-11-17, THE AUTHOR'S, docs/open_items/run_2026-11-17/
# DECISIONS.md; its details the supervisor's; record phase_P1_b1_ratio.md).  A
# default deuteron run now reads MillerB1() as the ratio b1/F1 at Q0^2 =
# B1_REFERENCE_Q2 = 2.5 GeV2: b1(x, Q2) = MillerB1(x) F1(x, Q2)/F1(x, Q0^2), the
# kernel's own F1, then the J1 note's confinement.  By convention 8 the headline is
# the b1 a default run carries, so since K1 `miller_shipped` is taken from the
# RUN-SURFACE kernel itself (`shipped_inclusive_kernel(deuteron(),
# b1_reference_q2=B1_REFERENCE_Q2)`, asserted equal to the default Pipeline's own
# kernel at the six bins in test_bench_spin_deuteron.py) at each bin's (<x>, <Q2>);
# below the grid floor its F1 is the frozen edge's (D7).  The row's rule, its
# tolerances and the other configurations are untouched (`miller_x1` stays twice
# the frozen table: registry row 2's normalisation option of the table, not a
# run).  The frozen-table reading shipped until K1 is kept as ONE RECORDED sub-row
# (FROZEN_SUBROW: chi2(b1) 5.262/6, chi2(A_zz) 7.308/6 on the shipped F1^d, and
# J4's in-grid numbers 1.894/3, 1.178/3 beside it); J4's in-grid sub-row is
# re-taken on the new headline.
#: the headline: the shipped b1 (A5, 2026-11-03; "cdks_conv_mstw" until then; K1 (2026-11-17): the run-surface kernel's b1)
HEADLINE = "miller_shipped"
#: K1 (2026-11-17): the recorded sub-row of the frozen table (no verdict)
FROZEN_SUBROW = "miller_shipped[frozen table: the reading shipped until K1] (K1, recorded)"
FROZEN_LABEL = ("the frozen MillerB1() table, b1 independent of Q2 -- the headline's "
                "reading until K1 (2026-11-17), kept for comparison: no verdict")
#: J4 (2026-11-10): the recorded in-grid sub-row's name (a diagnostic, no verdict)
INGRID_SUBROW = "miller_shipped[in-grid: <Q2> >= the shipped grid's Q2 floor] (J4, recorded)"
INGRID_LABEL = ("chosen after the values were seen (J4, 2026-11-10, the supervisor's "
                "choice under the author's delegation): a diagnostic, no verdict")


def measure():
    import lipolgen._lipolgen as _l
    rows = load_table()
    r1990 = r1990_factory()
    models, blocked, f2, f2_name = build_models()
    miller_undo = 1.0 / _l.B1_MILLER_TABLE_TO_PER_NUCLEON
    import lipolgen as lg  # noqa: WPS433
    shipped = lg.shipped_inclusive_kernel(lg.deuteron())   # A5: the shipped F2^d and R
    run_k = lg.shipped_inclusive_kernel(lg.deuteron(), b1_reference_q2=_l.B1_REFERENCE_Q2)  # K1 (2026-11-17): the run surface's b1
    # J4 (2026-11-10): the shipped grid's own Q2 floor, from the one shared object the kernel holds
    q2_floor = _l.make_lhapdf_unpol_sf(_l.UNPOL_SF_DEFAULT).q2_min
    per_bin = []
    for r in rows:
        x, q2 = r["x"], r["q2"]
        f2p, f2n = f2.f2p(x, q2), f2.f2n(x, q2)
        f2d = 0.5 * (f2p + f2n)                 # rule 2: the isoscalar nucleon
        g2 = 4.0 * M_NUCLEON ** 2 * x * x / q2
        r90, r98 = r1990(x, q2), _l.r1998(x, q2)
        f1d = (1.0 + g2) * f2d / (2.0 * x * (1.0 + r90))
        f1d_r98 = (1.0 + g2) * f2d / (2.0 * x * (1.0 + r98))
        t = shipped.tables(x, q2)                       # A5 (2026-11-03)
        r_sh = lg.r_sigma_lt(x, q2)
        f1d_sh = (1.0 + g2) * t.f2 / (2.0 * x * (1.0 + r_sh))
        b = {}
        for name, m in models.items():
            v = 0.0 if m is None else m.b1(x, q2, 0.0)  # f1 arg unused here
            b[name] = v * miller_undo if name == "miller_x1" else v
        frozen = b[HEADLINE]                            # K1 (2026-11-17): MillerB1(), the frozen table
        tk = run_k.tables(x, q2)
        b[HEADLINE] = tk.b1                             # K1: the run-surface kernel's b1 (scaled, confined)
        per_bin.append(dict(
            x=x, q2=q2, data=r, f2d=f2d, f1d=f1d, f1d_r1998=f1d_r98,
            r1990=r90, r1998=r98, f1_ratio_to_implied=f1d / r["f1_implied"],
            b1=b,
            azz={k: -(2.0 / 3.0) * v / f1d for k, v in b.items()},
            azz_r1998={k: -(2.0 / 3.0) * v / f1d_r98 for k, v in b.items()},
            f2d_shipped=t.f2, r_shipped=r_sh, f1d_shipped=f1d_sh,
            f1d_shipped_kernel=t.f1,
            azz_shipped={k: -(2.0 / 3.0) * v / f1d_sh for k, v in b.items()},
            q2_grid_floor=q2_floor, below_grid_floor=q2 < q2_floor,           # J4 (2026-11-10)
            f2d_frozen_at_floor=(q2 < q2_floor and t.f2 == shipped.tables(x, q2_floor).f2),
            f1d_shipped_over_implied=f1d_sh / r["f1_implied"],
            b1_frozen_table=frozen,                                            # K1 (2026-11-17)
            k1_factor=run_k.b1_reference_factor(x, tk.f1),                     # F1(x, <Q2>)/F1(x, 2.5)
            k1_f1_kernel=tk.f1, k1_unconfined_b1=run_k.unconfined_b1(x, q2, tk.f1),
            b1_frozen={HEADLINE: frozen},
            azz_frozen={HEADLINE: -(2.0 / 3.0) * frozen / f1d},
            azz_shipped_frozen={HEADLINE: -(2.0 / 3.0) * frozen / f1d_sh},
        ))
    return per_bin, blocked, f2_name


def ingrid_subrow(per_bin, name=HEADLINE, b1_key="b1", azz_key="azz_shipped"):
    """J4 (2026-11-10): the headline model on the bins at or above the shipped
    grid's Q2 floor -- RECORDED, chosen after the values were seen, no verdict.
    K1 (2026-11-17): `b1_key` / `azz_key` select the frozen table's numbers."""
    ing = [p for p in per_bin if not p["below_grid_floor"]]
    below = [p for p in per_bin if p["below_grid_floor"]]
    cb, wb = chi2_of(ing, name, "b1", key=b1_key)
    ca, wa = chi2_of(ing, name, "azz", key=azz_key)
    return dict(
        status="recorded", in_verdict=False, label=INGRID_LABEL,
        q2_grid_floor=per_bin[0]["q2_grid_floor"], unpol_sf="CT18NLO (UNPOL_SF_DEFAULT)",
        n_bins=len(ing), ndf=len(ing), chi2_b1=cb, within_b1=wb, chi2_azz=ca,
        within_azz=wa, azz_f1d="shipped, as the headline's",
        x_in_grid=[p["x"] for p in ing],
        bins_below_floor="; ".join("x = %g, <Q2> = %g GeV2" % (p["x"], p["q2"])
                                   for p in below),
        frozen_edge=("below the floor the shipped F2^d is the floor's (D7, frozen; "
                     "checked bin by bin: %s); shipped F1^d / the F1^d HERMES's own b1 "
                     "and A_zz columns imply = %s there, %s on the in-grid bins" % (
                         "all" if all(p["f2d_frozen_at_floor"] for p in below) else "NOT all",
                         " / ".join("%.3f" % p["f1d_shipped_over_implied"] for p in below),
                         " / ".join("%.3f" % p["f1d_shipped_over_implied"] for p in ing))))


def chi2_of(per_bin, name, obs, key=None):
    key = key or obs
    c, within = 0.0, 0
    for p in per_bin:
        d = p["data"][obs]
        e = p["data"][obs + "_err"]
        m = p[key][name]
        c += ((m - d) / e) ** 2
        within += abs(m - d) <= e
    return c, within


def run(verbose=True):
    per_bin, blocked, f2_name = measure()
    configs = {}
    for name in CONFIGS:
        if name in blocked:
            configs[name] = dict(status="blocked", reason=blocked[name])
            continue
        cb, wb = chi2_of(per_bin, name, "b1")
        ca, wa = chi2_of(per_bin, name, "azz")
        ca98, _ = chi2_of(per_bin, name, "azz", key="azz_r1998")
        extra = {}
        if name == HEADLINE:   # A5 (2026-11-03): A_zz on the shipped F1^d and R
            extra = dict(chi2_azz_mstw_r1990=ca, within_azz_mstw_r1990=wa,
                         azz_f1d="shipped: CT18NLO F2 (shipped deuteron kernel), "
                                 "r_sigma_lt, (1 + gamma^2)")
            ca, wa = chi2_of(per_bin, name, "azz", key="azz_shipped")
        configs[name] = dict(
            chi2_b1=cb, within_b1=wb, chi2_azz=ca, within_azz=wa,
            chi2_azz_r1998=ca98, ndf=len(per_bin),
            status="pass" if (wb == len(per_bin) and wa == len(per_bin)) else "fail",
            **extra)
    head = HEADLINE        # A5 (2026-11-03): "cdks_conv_mstw" (or cdks_conv_default) until then
    h = configs[head]
    d = configs["cdks_conv_mstw"]
    report = dict(
        name=NAME,
        generator="%s (the shipped b1; the headline under A5; A_zz on the shipped "
                  "F1^d and R): chi2/ndf b1 = %.2f/%d, A_zz = %.2f/%d; %d/%d b1 bins "
                  "within tol; cdks_conv_mstw (declared until A5, recorded): %s" % (
                      head, h["chi2_b1"], h["ndf"], h["chi2_azz"],
                      h["ndf"], h["within_b1"], h["ndf"],
                      "blocked" if d.get("status") == "blocked" else
                      "b1 %.2f/%d, A_zz %.2f/%d" % (d["chi2_b1"], d["ndf"],
                                                     d["chi2_azz"], d["ndf"])),
        reference="HERMES PRL 95 (2005) 242001 Table II, 6 bins, per nucleon",
        tolerance="per bin sqrt(stat^2+syst^2)",
        status=h["status"], headline=head, configs=configs,
        per_bin=per_bin, f2_source=f2_name,
        subrows={INGRID_SUBROW: ingrid_subrow(per_bin),      # J4 (2026-11-10): recorded
                 FROZEN_SUBROW: frozen_subrow(per_bin)},     # K1 (2026-11-17): recorded
    )
    if verbose:
        print("F1^d rebuilt as HERMES did, with F2 = (F2p+F2n)/2 from %s and "
              "R1990, vs the F1^d implied by their table (-2/3 b1/A_zz):" % f2_name)
        for p in per_bin:
            print("  x=%.3f Q2=%.2f  F1_rebuilt=%.4f F1_implied=%.4f ratio=%.4f"
                  "  R1990=%.3f R1998=%.3f (1+R90)/(1+R98)=%.4f" % (
                      p["x"], p["q2"], p["f1d"], p["data"]["f1_implied"],
                      p["f1_ratio_to_implied"], p["r1990"], p["r1998"],
                      (1 + p["r1990"]) / (1 + p["r1998"])))
        print("b1^d per nucleon, units 1e-2 (HERMES +- quadrature error):")
        print("  x      HERMES          " + " ".join("%17s" % n for n in CONFIGS
                                                     if n not in blocked))
        for p in per_bin:
            d = p["data"]
            print("  %.3f  %+6.2f +- %5.2f  " % (p["x"], 100 * d["b1"], 100 * d["b1_err"])
                  + " ".join("%17.4f" % (100 * p["b1"][n]) for n in CONFIGS
                             if n not in blocked))
        print("A_zz^d, units 1e-2:")
        for p in per_bin:
            d = p["data"]
            print("  %.3f  %+6.2f +- %5.2f  " % (p["x"], 100 * d["azz"], 100 * d["azz_err"])
                  + " ".join("%17.4f" % (100 * p["azz"][n]) for n in CONFIGS
                             if n not in blocked))
        for n in CONFIGS:
            c = configs[n]
            if c["status"] == "blocked":
                print("REPORT | %s[%s] | BLOCKED: %s" % (NAME, n, c["reason"]))
                continue
            print("REPORT | %s[%s] | chi2/ndf b1 = %.3f/%d (%d/6 within), A_zz = "
                  "%.3f/%d (%d/6 within; with R1998 %.3f) | HERMES Table II | "
                  "tol per-bin quadrature | %s" % (
                      NAME, n, c["chi2_b1"], c["ndf"], c["within_b1"],
                      c["chi2_azz"], c["ndf"], c["within_azz"],
                      c["chi2_azz_r1998"], c["status"].upper()))
        g = report["subrows"][INGRID_SUBROW]                  # J4 (2026-11-10)
        print("J4 in-grid sub-row, %s, %s: <Q2> >= %.6g GeV2 (%d bins): chi2/ndf b1 = "
              "%.3f/%d (%d within), A_zz = %.3f/%d (%d within, shipped F1^d); below the "
              "floor: %s; %s" % (HEADLINE, g["label"], g["q2_grid_floor"], g["n_bins"],
                                 g["chi2_b1"], g["ndf"], g["within_b1"], g["chi2_azz"],
                                 g["ndf"], g["within_azz"], g["bins_below_floor"],
                                 g["frozen_edge"]))
        f = report["subrows"][FROZEN_SUBROW]                  # K1 (2026-11-17)
        print("K1 frozen-table sub-row, %s: chi2/ndf b1 = %.3f/%d (%d within), A_zz = "
              "%.3f/%d (%d within, shipped F1^d; %.3f on MSTW/R1990); in-grid b1 %.3f/%d, "
              "A_zz %.3f/%d" % (f["label"], f["chi2_b1"], f["ndf"], f["within_b1"],
                                f["chi2_azz"], f["ndf"], f["within_azz"],
                                f["chi2_azz_mstw_r1990"], f["in_grid"]["chi2_b1"],
                                f["in_grid"]["ndf"], f["in_grid"]["chi2_azz"],
                                f["in_grid"]["ndf"]))
        for p in per_bin:
            print("  K1 x=%.3f <Q2>=%.2f  F1(x,<Q2>)/F1(x,2.5) = %.6f  b1 frozen %.6e -> "
                  "headline %.6e" % (p["x"], p["q2"], p["k1_factor"], p["b1_frozen_table"],
                                     p["b1"][HEADLINE]))
    return report


def frozen_subrow(per_bin, name=HEADLINE):
    """K1 (2026-11-17): the frozen MillerB1() table -- the headline's reading
    until K1 -- on all six bins, A_zz on the shipped F1^d as the headline's,
    with J4's in-grid numbers on it beside; RECORDED, no verdict."""
    cb, wb = chi2_of(per_bin, name, "b1", key="b1_frozen")
    ca, wa = chi2_of(per_bin, name, "azz", key="azz_shipped_frozen")
    cm, wm = chi2_of(per_bin, name, "azz", key="azz_frozen")
    g = ingrid_subrow(per_bin, name, b1_key="b1_frozen", azz_key="azz_shipped_frozen")
    return dict(
        status="recorded", in_verdict=False, label=FROZEN_LABEL, ndf=len(per_bin),
        chi2_b1=cb, within_b1=wb, chi2_azz=ca, within_azz=wa,
        chi2_azz_mstw_r1990=cm, within_azz_mstw_r1990=wm,
        azz_f1d="shipped, as the headline's",
        in_grid=dict(n_bins=g["n_bins"], ndf=g["ndf"], chi2_b1=g["chi2_b1"],
                     within_b1=g["within_b1"], chi2_azz=g["chi2_azz"],
                     within_azz=g["within_azz"], label=g["label"]))


# ------------------------------------------------------------------ figure

#: the figure's palette slot, marker and legend name of each configuration
#: (the null baseline b1 = 0 is drawn as the zero line, not as a series)
PLOT_SLOTS = {
    "miller_shipped": ("default", "D", "Miller x 0.5 (shipped; b1/F1 at 2.5 GeV2, K1)"),
    "miller_x1": ("alt3", "o", "Miller x 1 (row 2's option)"),
    "cdks_conv_mstw": ("mstw", "s", "CDKS conv., MSTW (A = 2 gate)"),
    "cdks_conv_default": ("alt", "^", "CDKS conv., ToyF2"),
    "cdks_digitized": ("alt2", "v", "CDKS Fig. 4, digitized"),
}


def plot(report, out_dir=None):
    """x b1^d per nucleon in HERMES's six bins (ink, stat (+) syst), each
    configuration's b1 at the bin's own (<x>, <Q2>) as a marker in its slot
    colour (side by side), the null baseline b1 = 0 as the zero line, over
    the per-bin pulls (model - data)/sigma with the tolerance |pull| <= 1
    shaded.  x b1 rather than b1 so that all six bins are readable; x is
    the bin's, so each pull is b1's.  A_zz is in the row, not drawn."""
    if HERE not in sys.path:
        sys.path.insert(0, HERE)
    from _plots import BAND, MUTED, new_figure, save, status_title, style_series
    bins, configs = report["per_bin"], report["configs"]
    xs = [p["x"] for p in bins]
    fig, (ax, axp) = new_figure(two_panel=True)
    axp.axhspan(-1.0, 1.0, color=BAND, alpha=0.8, linewidth=0, zorder=1)
    ax.axhline(0.0, color=MUTED, linewidth=1.0, linestyle="--", zorder=2,
               label="b1 = 0 (null baseline): %.2f/%d"
               % (configs["null_b1_zero"]["chi2_b1"], configs["null_b1_zero"]["ndf"]))
    ax.errorbar(xs, [p["x"] * p["data"]["b1"] for p in bins],
                yerr=[p["x"] * p["data"]["b1_err"] for p in bins],
                **style_series("reference", label="HERMES 2005 Table II (stat (+) syst)"))
    lo_q2 = [p for p in bins if p.get("below_grid_floor")]   # J4 (2026-11-10): the frozen edge
    if lo_q2:
        ax.plot([p["x"] for p in lo_q2], [p["x"] * p["data"]["b1"] for p in lo_q2], "s",
                markersize=11, markerfacecolor="none", markeredgecolor=MUTED, zorder=6,
                label="<Q2> below the shipped grid's floor %.3g GeV2 (F2 frozen, D7; J4)"
                % lo_q2[0]["q2_grid_floor"])
    for k, (n, (slot, mk, label)) in enumerate(PLOT_SLOTS.items()):
        c = configs[n]
        if c["status"] == "blocked":
            continue
        dodge = 10 ** (0.03 * (k - 2))        # side by side at each bin
        ax.plot([x * dodge for x in xs], [p["x"] * p["b1"][n] for p in bins],
                **style_series(slot, label="%s: %.2f/%d, %d/%d in"
                               % (label, c["chi2_b1"], c["ndf"], c["within_b1"],
                                  c["ndf"]), marker=mk))
        pk = style_series(slot, marker=mk)
        pk.update(label="_nolegend_", markersize=5)
        axp.plot([x * dodge for x in xs],
                 [(p["b1"][n] - p["data"]["b1"]) / p["data"]["b1_err"] for p in bins],
                 **pk)
    fz = report.get("subrows", {}).get(FROZEN_SUBROW)      # K1 (2026-11-17): the frozen table
    if fz:
        ax.plot(xs, [p["x"] * p["b1_frozen_table"] for p in bins], "x", color=MUTED,
                markersize=7, zorder=5, label="Miller x 0.5, frozen table (until K1): "
                "%.2f/%d, %d/%d in" % (fz["chi2_b1"], fz["ndf"], fz["within_b1"], fz["ndf"]))
    ax.set_xscale("log")
    ax.set_ylabel("x b1^d (per nucleon)")
    lo, hi = ax.get_ylim()
    ax.set_ylim(lo - 0.75 * (hi - lo), hi)  # room for the legend below the data
    ax.legend(loc="lower center", ncol=2, fontsize=7,
              title="chi2(b1)/ndf, bins within tolerance", title_fontsize=7)
    axp.set_xlabel("<x> of the HERMES bin (each at its own <Q2>, 0.51-4.69 GeV^2)")
    axp.set_ylabel("(model - data)/sigma")
    h = configs[report["headline"]]
    status_title(fig, NAME, report["status"],
                 "tolerance: every bin within stat (+) syst (grey band); headline "
                 "%s %.2f/%d, %d/%d bins within" % (
                     report["headline"], h["chi2_b1"], h["ndf"], h["within_b1"],
                     h["ndf"]))
    return save(fig, NAME, out_dir)


if __name__ == "__main__":
    # Exit convention (module docstring, validation/benchmarks/README.md):
    # 0 = ran and printed its REPORT row, whatever the verdict; 2 = broken.
    # `--plot` also writes the row's figure (plot()); the same rule holds.
    try:
        rep = run()
        if "--plot" in sys.argv[1:]:
            print("figure: %s" % plot(rep))
    except Exception:  # any exception means the harness itself is broken
        traceback.print_exc()
        sys.exit(2)
    sys.exit(0)
