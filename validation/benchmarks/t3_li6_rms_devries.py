#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Scout item 6 / 03_data_nuclear.md sec. 2 ("The de Vries 6Li rms is a real,
immediate gate") -- the 6Li charge rms radius of the tree's C0 form factor
against the electron-scattering analyses (tier T3).

WHAT IS COMPARED.  The charge rms radius of the SHIPPED
`HoSpin1FF.for_ion(li6())` (C0 shape `ho`, `fold_nucleon = True`, nothing
moved), read off its q -> 0 slope,

    <r^2> = -6 d[F_c(q)/F_c(0)]/d(q^2) at q -> 0,      q^2 = t/(hbar c)^2,

against the 6Li rows of de Vries, de Jager, de Vries, ADNDT 36 (1987) 495,
Table I -- 2.54(5), 2.56(5), 2.57(10) fm, three electron-scattering
analyses (Su67, Li71a, Bu72) -- and Angeli & Marinova, ADNDT 99 (2013) 69,
2.589(39) fm.  All four are vendored in `data/devries1987_li6_rms.json` as
a SURVEY TRANSCRIPTION (03_data_nuclear.md lines 182-188 and 207; rc.hpp
line 543): the ADNDT pages were not re-read here.

NUCLEON FINITE SIZE IS FOLDED IN, like for like with a charge radius.  The
shipped `fc` is F_c(t) = Z F_point(q) [G_E^p(t) + G_E^n(t)] (rc.hpp,
`HoSpin1FFOptions::fold_nucleon = true`): the ISOSCALAR sum (N = Z = 3),
G_E^p the dipole 1/(1 + t/0.71)^2 and G_E^n Galster's, both from rc.hpp's
`nucleon_ff`.  By the product rule its slope radius is EXACTLY
    r^2 = r^2_point + r^2_p + (N/Z) r^2_n,     N/Z = 1,
with r^2_p = 12/0.71 GeV^-2 (hbar c)^2 = 0.6581 fm^2 (the dipole's, not the
0.7071 of rc.hpp's own recipe) and r^2_n = -0.1269 fm^2 (Galster's slope,
not the recipe's -0.1155).  NO Darwin-Foldy term and no spin-orbit term are
in `fc`.  Each piece is measured separately below (the point shape with
`fold_nucleon = false`, the nucleon slopes from `lg.nucleon_ff`) and the
decomposition is asserted by the pytest.

THE SLOPE, NUMERICALLY.  g(s) = 6 [1 - F(s)/F(0)]/s at s = q^2 = S0/2^i,
i = 0..4 (S0 = 0.05 fm^-2), Richardson-extrapolated to s -> 0 (g is analytic
in s: g = <r^2> - c_1 s + c_2 s^2 - ...).  The numerical error carried in
the row is |T_44 - T_33| of the Richardson table (7.4e-9 fm^2, i.e.
1.4e-9 fm in r -- a conservative bound); round-off in g is ~2e-13 fm^2 at
the smallest s, a few 1e-12 after extrapolation.  Cross-checks: the HO
closed form (3/2) a^2 (2 + 5 alpha)/(2 + 3 alpha) + r^2_p + r^2_n (measured
residual -3.8e-12 fm^2), and the one-sided eps = 1e-3 fm^-1 difference
t3_li6_charge_ff_fb.py prints (2.5710001 fm, 7e-7 fm low).

TOLERANCE, DECLARED BEFORE THE VERDICT -- and AFTER the number was known:
the scout's preview printed r_ch = 2.5710 fm, and t3_li6_charge_ff_fb.py
already prints the same number, before this tolerance was chosen.  It is:

    PASS iff |r_tree - r_i| <= 2 sigma_i for EACH of the three de Vries
    Table I analyses (the electron-scattering determinations).

Angeli-Marinova is compared and its pull printed, but it is NOT in the
verdict: `LI6_R2_POINT_FM2` = 6.0788 fm^2 was DERIVED from its 2.589 fm
(rc.hpp: minus <r^2>_p = 0.7071, plus -<r^2>_n = 0.1155, minus the
Darwin-Foldy 0.033), so r_tree - 2.589 measures the tree's own folding
bookkeeping, not an independent datum.  The row prints that closure,
decomposed: (dipole r^2_p - 0.7071) + (Galster r^2_n + 0.1155) - 0.033
(the Darwin-Foldy term `fc` does not add back) + ((a, alpha) point r^2 -
the recipe's) = r^2_tree - 2.589^2 exactly.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations -- NOT APPLICABLE: no b1, no A-scaling assumption.
  2. isoscalar denominator -- NOT APPLICABLE as a structure-function
     denominator (none enters).  Its analogue here is the nucleon folding,
     which is ISOSCALAR, G_E^p + G_E^n with N = Z: asserted by the pytest
     through the decomposition above (a proton-only fold would give
     2.5956 fm, not 2.5710).
  3. per nucleus / per nucleon -- both sides PER NUCLEUS: the de Vries and
     Angeli radii are of the charge distribution normalised to Z = 3; the
     tree's F_c(0) = Z = 3.  The slope is taken of F_c/F_c(0), so Z cancels
     (asserted: F_c(0) = 3, and a rescaled F gives the same radius).

WHAT THIS DOES NOT VALIDATE.  The SHAPE of F_C beyond its slope (its first
zero is row 5, t3_li6_charge_ff_fb.py, a recorded FAIL at +0.27 fm^-1); the
C0 band (`C0Shape.VmcFt` carries the same <r^2>_point by construction, so
the radius cannot tell the two edges apart -- printed, not compared); the
magnetization radius (row 6, BLOCKED) or the quadrupole moment (an input);
7Li (`HoSpin1FF.for_ion` throws for it); and it is NOT a prediction of a
nuclear model: the point radius is an INPUT of `HoSpin1FF`, derived from
Angeli.  A PASS says that input, folded as the tree folds it, sits inside
the independent electron-scattering determinations -- no more.

Run:  source env.sh && python3 validation/benchmarks/t3_li6_rms_devries.py

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
import time
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data", "devries1987_li6_rms.json")
#: sha256 of every byte of the vendored file, and of its canonical body
#: (json.dumps(body, sort_keys=True, separators=(",", ":"))), which the file
#: also records.  Both are re-checked on EVERY read.
FILE_SHA256 = "ffab49a1b67d92a209de3341de994c6b1fc115fe7f69044ef054d3154a53ad3c"
BODY_SHA256 = "c125dd346bfd9ea6e39f65b7b8ac8dccaae8ba20eb9bf1ad6a2052318b7238cf"

NAME = "t3_li6_rms_devries"
MAX_PULL = 2.0          # tolerance: |r - r_i| <= 2 sigma_i, each de Vries row
S0_FM2 = 0.05           # Richardson: s = q^2 = S0/2^i fm^-2, i = 0..LEVELS-1
LEVELS = 5
EPS_ONE_SIDED = 1e-3    # fm^-1, t3_li6_charge_ff_fb.py's rms_folded step
#: rc.hpp:541-546's recipe for LI6_R2_POINT_FM2, as typed in its prose
#: [fm^2] -- used ONLY to decompose the Angeli closure, never as an input.
RECIPE_R2_P = 0.7071
RECIPE_R2_N = -0.1155
RECIPE_DARWIN_FOLDY = 0.033


# ---------------------------------------------------------------- reference

def read_vendored(path=DATA):
    """The vendored document -- RuntimeError unless the file's bytes AND its
    canonical body hash to the recorded digests."""
    with open(path, "rb") as f:
        raw = f.read()
    got = hashlib.sha256(raw).hexdigest()
    if got != FILE_SHA256:
        raise RuntimeError("%s: vendored file sha256 %s != recorded %s"
                           % (path, got, FILE_SHA256))
    doc = json.loads(raw.decode("ascii"))
    canon = json.dumps(doc["body"], sort_keys=True, separators=(",", ":"),
                       ensure_ascii=True, allow_nan=False)
    body = hashlib.sha256(canon.encode("ascii")).hexdigest()
    if not body == doc["_body_sha256"] == BODY_SHA256:
        raise RuntimeError("%s: body sha256 %s; file records %s, harness %s"
                           % (path, body, doc["_body_sha256"], BODY_SHA256))
    return doc


def load_reference(path=DATA):
    """(the three de Vries Table I rows, the Angeli-Marinova point)."""
    pts = [dict(p, rms_fm=float(p["rms_fm"]), sigma_fm=float(p["sigma_fm"]))
           for p in read_vendored(path)["body"]["points"]]
    dv = [p for p in pts if p["compilation"].startswith("de Vries")]
    am = [p for p in pts if p["compilation"].startswith("Angeli")]
    if len(pts) != 4 or len(dv) != 3 or len(am) != 1:
        raise RuntimeError("%s: expected 3 de Vries rows and 1 Angeli-Marinova "
                           "point, found %d/%d of %d" % (path, len(dv), len(am),
                                                          len(pts)))
    return dv, am[0]


# ---------------------------------------------------------------- the slope

def richardson_r2(f, s0=S0_FM2, levels=LEVELS, f0=None):
    """<r^2> = -6 dF/ds / F(0) at s -> 0 for F(s), s = q^2 [fm^-2].
    Returns (<r^2>, |T_nn - T_(n-1)(n-1)|, the diagonal)."""
    f0 = f(0.0) if f0 is None else f0

    def g(s):
        return 6.0 * (1.0 - f(s) / f0) / s

    tab = [[g(s0 / 2 ** i)] for i in range(levels)]
    for i in range(1, levels):
        for j in range(1, i + 1):
            tab[i].append(tab[i][j - 1]
                          + (tab[i][j - 1] - tab[i - 1][j - 1]) / (2 ** j - 1))
    diag = [tab[i][i] for i in range(levels)]
    return diag[-1], abs(diag[-1] - diag[-2]), diag


def slope_of_zero_at_origin(f, s0=S0_FM2, levels=LEVELS):
    """-6 dF/ds at s -> 0 for an F with F(0) = 0 (G_E^n), same scheme."""
    return richardson_r2(lambda s: 1.0 + f(s), s0, levels, f0=1.0)


def tree_side():
    import lipolgen as lg  # noqa: WPS433 (env.sh puts it on the path)
    hc2 = lg.HBARC_GEV_FM ** 2
    ion = lg.li6()
    ff = lg.HoSpin1FF.for_ion(ion)                    # SHIPPED, nothing moved
    opt = ff.options
    if not opt.fold_nucleon or opt.c0_shape != lg.C0Shape.Ho:
        raise RuntimeError("HoSpin1FF.for_ion(li6()) is no longer the shipped "
                           "(ho, fold_nucleon) configuration this row names")
    fc = lambda s: ff.fc(s * hc2)                     # noqa: E731
    r2, err, diag = richardson_r2(fc)
    fc0 = fc(0.0)
    eps2 = EPS_ONE_SIDED ** 2
    r2_one_sided = -6.0 * (fc(eps2) / fc0 - 1.0) / eps2

    o_pt = lg.HoSpin1FFOptions()
    o_pt.fold_nucleon = False
    ff_pt = lg.HoSpin1FF.for_ion(ion, o_pt)
    r2_pt, err_pt, _ = richardson_r2(lambda s: ff_pt.fc(s * hc2))
    r2_p, err_p, _ = richardson_r2(lambda s: lg.nucleon_ff(s * hc2).ge_p)
    r2_n, err_n, _ = slope_of_zero_at_origin(lambda s: lg.nucleon_ff(s * hc2).ge_n)
    a, alpha = opt.a_fm, opt.alpha
    r2_pt_closed = 1.5 * a * a * (2 + 5 * alpha) / (2 + 3 * alpha)

    o_vmc = lg.HoSpin1FFOptions()
    o_vmc.c0_shape = lg.C0Shape.VmcFt
    ff_vmc = lg.HoSpin1FF.for_ion(ion, o_vmc)
    r2_vmc, err_vmc, _ = richardson_r2(lambda s: ff_vmc.fc(s * hc2))

    r = math.sqrt(r2)
    return dict(
        r2=r2, r=r, r2_error=err, r_error=err / (2.0 * r), richardson_diag=diag,
        r_one_sided=math.sqrt(r2_one_sided), fc0=fc0, z=ion.Z,
        r2_point=r2_pt, r2_point_error=err_pt, r2_point_closed=r2_pt_closed,
        r2_point_constant=lg.LI6_R2_POINT_FM2,
        r2_p=r2_p, r2_p_error=err_p, r2_n=r2_n, r2_n_error=err_n,
        decomposition_residual=r2 - (r2_pt + r2_p + r2_n),
        closed_form_residual=r2 - (r2_pt_closed + r2_p + r2_n),
        r_proton_only_fold=math.sqrt(r2_pt + r2_p),
        r_vmcft_edge=math.sqrt(r2_vmc), r2_vmcft_error=err_vmc,
        a_fm=a, alpha=alpha, provenance=ff.provenance(),
        config=dict(c0_shape="ho", fold_nucleon=True, s0_fm2=S0_FM2,
                    levels=LEVELS),
    )


# ---------------------------------------------------------------- verdict

def verdict(r, de_vries):
    """Per de Vries row: pull and |pull| <= MAX_PULL.  PASS iff all three."""
    rows = []
    for p in de_vries:
        pull = (r - p["rms_fm"]) / p["sigma_fm"]
        rows.append(dict(p, pull=pull, within=abs(pull) <= MAX_PULL,
                         within_1sigma=abs(pull) <= 1.0))
    return rows, ("pass" if all(x["within"] for x in rows) else "fail")


def angeli_closure(tree, angeli):
    """r^2_tree - r^2_Angeli, decomposed into what the tree's folding does
    differently from rc.hpp's recipe.  Recorded, never in the verdict."""
    r2_a = angeli["rms_fm"] ** 2
    recipe_pt = r2_a - RECIPE_R2_P - RECIPE_R2_N - RECIPE_DARWIN_FOLDY
    parts = dict(
        proton=tree["r2_p"] - RECIPE_R2_P,
        neutron=tree["r2_n"] - RECIPE_R2_N,
        darwin_foldy=-RECIPE_DARWIN_FOLDY,
        point=tree["r2_point"] - recipe_pt)
    return dict(
        pull=(tree["r"] - angeli["rms_fm"]) / angeli["sigma_fm"],
        distance_fm=tree["r"] - angeli["rms_fm"],
        dr2=tree["r2"] - r2_a, parts=parts,
        parts_residual=tree["r2"] - r2_a - sum(parts.values()),
        recipe_r2_point=recipe_pt)


def measure(path=DATA):
    de_vries, angeli = load_reference(path)
    tree = tree_side()
    rows, status = verdict(tree["r"], de_vries)
    lo = max(p["rms_fm"] - MAX_PULL * p["sigma_fm"] for p in de_vries)
    hi = min(p["rms_fm"] + MAX_PULL * p["sigma_fm"] for p in de_vries)
    return dict(tree=tree, de_vries=rows, angeli=angeli,
                closure=angeli_closure(tree, angeli), status=status,
                window=(lo, hi))


def run(path=DATA, verbose=True):
    t0 = time.perf_counter()
    m = measure(path)
    t, c = m["tree"], m["closure"]
    generator = ("HoSpin1FF.for_ion(li6()) C0 slope, nucleon-folded (G_E^p + "
                 "G_E^n): r_ch = %.5f fm (+- %.1e numerical)" % (t["r"], t["r_error"]))
    reference = ("de Vries ADNDT 36 (1987) Table I %s fm; Angeli-Marinova %s fm "
                 "(closure, recorded); survey transcription" % (
                     " / ".join(p["as_printed"] for p in m["de_vries"]),
                     m["angeli"]["as_printed"]))
    tolerance = ("abs(r - r_i) <= %g sigma_i for each de Vries row (i.e. r in "
                 "[%.2f, %.2f] fm)" % (MAX_PULL, m["window"][0], m["window"][1]))
    import lipolgen as lg  # noqa: WPS433
    report = dict(m, name=NAME, generator=generator, reference=reference,
                  tolerance=tolerance,
                  config=dict(lipolgen=getattr(lg, "__version__", None),
                              python=sys.version.split()[0],
                              data_file_sha256=FILE_SHA256))
    report["runtime_s"] = time.perf_counter() - t0
    if verbose:
        print("tree: F_c(0) = %.12g (Z = %d); Richardson diagonal <r^2> = %s fm^2"
              % (t["fc0"], t["z"], ", ".join("%.10f" % v for v in t["richardson_diag"])))
        print("tree: r_ch = %.9f fm, <r^2> = %.10f fm^2 (+- %.1e); one-sided "
              "eps = %g fm^-1: %.7f fm" % (t["r"], t["r2"], t["r2_error"],
                                           EPS_ONE_SIDED, t["r_one_sided"]))
        print("decomposition: <r^2>_point %.8f (closed form (a, alpha) = (%.4f, "
              "%.5f): %.8f; LI6_R2_POINT_FM2 %.4f) + r^2_p %.8f (dipole) + r^2_n "
              "%.8f (Galster) = %.8f; residual %.1e, vs closed form %.1e"
              % (t["r2_point"], t["a_fm"], t["alpha"], t["r2_point_closed"],
                 t["r2_point_constant"], t["r2_p"], t["r2_n"],
                 t["r2_point"] + t["r2_p"] + t["r2_n"],
                 t["decomposition_residual"], t["closed_form_residual"]))
        print("recorded: point radius %.6f fm; proton-only fold would give "
              "%.6f fm; the C0 band's VmcFt edge %.6f fm (same <r^2>_point by "
              "construction)" % (math.sqrt(t["r2_point"]), t["r_proton_only_fold"],
                                 t["r_vmcft_edge"]))
        for p in m["de_vries"]:
            print("  de Vries %-5s %-4s %s fm (q %s fm^-1): pull %+.2f -> %s"
                  % (p["source_key"], p["model"], p["as_printed"],
                     "-".join("%.2f" % v for v in p["q_range_fm-1"]), p["pull"],
                     "within %g sigma" % MAX_PULL if p["within"] else "OUTSIDE"))
        pp = c["parts"]
        print("  Angeli-Marinova %s fm: pull %+.2f (distance %+.4f fm) -- CLOSURE, "
              "not in the verdict: r^2_tree - 2.589^2 = %+.5f fm^2 = %+.5f "
              "(dipole r^2_p - 0.7071) %+.5f (Galster r^2_n + 0.1155) %+.5f "
              "(no Darwin-Foldy in fc) %+.5f ((a, alpha) point r^2 - recipe "
              "%.6f); residual %.1e" % (
                  m["angeli"]["as_printed"], c["pull"], c["distance_fm"], c["dr2"],
                  pp["proton"], pp["neutron"], pp["darwin_foldy"], pp["point"],
                  c["recipe_r2_point"], c["parts_residual"]))
        for p in m["de_vries"]:
            print("SUBROW | %s[%s] | r_ch = %.5f fm | de Vries Table I %s %s fm | "
                  "abs(pull) <= %g | %s (pull %+.2f)"
                  % (NAME, p["source_key"], t["r"], p["model"], p["as_printed"],
                     MAX_PULL, "PASS" if p["within"] else "FAIL", p["pull"]))
        print("SUBROW | %s[angeli_closure] | r_ch = %.5f fm | Angeli-Marinova %s fm "
              "| not asserted (the tree's input) | RECORDED (pull %+.2f)"
              % (NAME, t["r"], m["angeli"]["as_printed"], c["pull"]))
        print("REPORT | %s | %s | %s | %s | %s" % (
            NAME, generator, reference, tolerance, m["status"].upper()))
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
