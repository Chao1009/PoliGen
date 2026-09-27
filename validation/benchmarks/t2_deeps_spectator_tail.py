#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""NEW row (proposed for BENCHMARK_PLAN.md sec. 4, not one of the ten) --
CLAS "Deeps" D(e,e'p_s)X: the SHAPE of the high-momentum spectator tail,
p_s = 0.30-0.53 GeV/c, against the tree's deuteron momentum distribution
(tier T2).

THE DATA.  Klimenko et al. (CLAS), PRC 73 (2006) 035212 (nucl-ex/0510032),
F2N(x*, Q2) x P(p_s, cos theta_pq), per deuteron, as the CLAS Physics
Database lists it (eid 90), vendored in `data/clas_deeps_klimenko2006_f2P.json`
(sha256 recorded below -- the SAME digest t2_bonus_spectator_shape.py
records -- and re-checked on every read).  The census of that file is
t2_bonus_spectator_shape.py's docstring; what THIS row does with it:
  * it uses ONLY the 60 blocks at the paper's Q2 bins 1.8 and 2.8 GeV^2
    (2 Q2 x 6 W* x 5 p_s, one block each -- asserted unique).  The 55 blocks
    the database labels Q2 = 0.18-0.65 (unverified labels, 25 duplicated
    label triples, one empty block) are NOT read by any comparison;
  * rows with value = stat = 0 are MISSING (dropped, 37 of them in the 60
    blocks), not measured zeros.  The only 2 negative values (mid 20 and 31,
    cos theta_pq = -0.85) also carry stat = 0 -- no statistical error, so
    no measurement a chi2 can weigh -- and are dropped with them: every row
    with stat = 0 is dropped (39), and none with stat > 0;
  * p_s is the paper's bin AVERAGE (0.30, 0.34, 0.39, 0.46, 0.53 GeV/c), not
    the database's label (its last bin says 560 MeV);
  * errors: stat (+) syst in quadrature.  The database gives no split of
    syst into point-to-point and normalisation parts; a normalisation part
    cancels in a shape and is here counted as if it did not -- conservative
    for p (it can only raise it).  The stat-only chi2 is recorded beside it.

WHY A SHAPE, AND WHY F2N CANCELS (sec. 8 rule 3).  At fixed (Q2, W*) the
struck neutron's x* is fixed -- up to the off-shellness p*^2 of the struck
nucleon, which moves x* by a few per cent across 0.30-0.53 GeV/c at fixed
W* and is NOT corrected (the tree has no F2n model at these W*, and W* =
0.94 is the quasi-elastic peak where "F2N" is the elastic term).  So within
a group of FIXED (Q2, W*, cos theta_pq), the p_s dependence of F2N x P is
the p_s dependence of P alone.  Each group gets ONE free scale, fitted
(`scaled_chi2`: s = sum d m / s^2 / sum m^2 / s^2, ndf = n - 1), which
absorbs F2N, the per-deuteron normalisation of P and the tree's n(k)
normalisation together; nothing absolute is compared.  A group with fewer
than 2 points carries no shape and is skipped (counted).

THE PWIA WINDOW -- the verdict.  cos theta_pq <= -0.3 (bin centres -0.85 ..
-0.35), where the paper finds PWIA adequate.  The model is the tree's
deuteron momentum density n(k) (unpolarized: the mean over the ion's
M = +1, 0, -1 of `TaggedModel(deuteron_channel(...)).n_of_kc(M, k, c)`,
which is isotropic, checked), taken to the data's LIGHT-CONE P by the
Frankfurt-Strikman prescription CODED IN THIS HARNESS (it is not in the
tree, and neither is the paper's own Eq. (8) with its flux factor):
    alpha_s = (E_s - p_s cos theta_pq) / (M_d / 2),  p_T = p_s sin theta_pq,
    k^2 = (m^2 + p_T^2) / (alpha_s (2 - alpha_s)) - m^2,
    P_LC(p_s, cos) = n(k) E_k / (2 - alpha_s),      E_k = sqrt(m^2 + k^2),
with E_s the on-shell proton energy, m = M_NUCLEON, M_d = nuclear_mass(1, 2).
Backward (alpha_s > 1) this k is much LARGER than p_s (k = 0.37 at p_s =
0.30, 0.95 GeV at p_s = 0.53, cos = -1), so the prescription is not a
detail.  RECORDED beside it, never in the verdict: the INSTANT-FORM reading
P = n(k = p_s), i.e. what the generator itself does (the tree puts the
spectator's rest-frame momentum AT k, `FsiKinematics` doc), and the LC model
times the tree's Glauber FSI weight.

Waves: the shipped default, `deuteron_channel()` = the analytic Hulthen
pair at P_D = 0.045 (the tagged-channel default) -- THE HEADLINE VERDICT;
AV18 (`source = ClusterWaveSource.VmcAV18`, fdeut.av18) and CD-Bonn
(`cdbonn_wave()`, Machleidt 2001 App. D, n ~ psi_s^2 + psi_d^2; the tagged
channel has no CD-Bonn source, so this sub-row reads the wave directly) as
sub-rows, judged by the same tolerance but NOT in the verdict.

TOLERANCE, declared BEFORE any chi2 of this row was computed (2026-09-27):
p(chi2, ndf) >= 0.01 on the COMBINED PWIA-window shape chi2 of the default
wave with the LC prescription.  The per-(Q2, W*) and per-Q2 chi2 are
reported, not asserted.

THE TRANSVERSE WINDOW -- recorded, NEVER asserted.  |cos theta_pq| < 0.3
(bin centres -0.25 .. 0.25), where the paper finds FSI dominant.  Per
(Q2, W*) the LC model is normalised by ONE scale fitted over that group's
whole PWIA window; data / PWIA at the transverse points is compared with
the tree's direction, `GlauberFsiWeight(model).weight(FsiKinematics(k = p_s,
cos_theta_k = cos theta_pq, w = W*, q2, x = x*))` (default options:
sigma_XN = 40 mb, eps = -0.5, ramp off; spectator (Z, A) = (1, 1), channel
TaggedDeuteronP), rescaled by the same backward-window fit of PWIA x weight
so both sides carry the same normalisation.  Per p_s bin: the weighted mean
data/PWIA, the tree's, and whether both sit on the same side of 1.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations -- NOT APPLICABLE (unpolarized).
  2. isoscalar denominator -- NOT APPLICABLE; its analogue, F2 of the struck
     NEUTRON, is held fixed (fixed Q2, W*) and absorbed in a free scale.
  3. per nucleus / per nucleon -- the data are a reduced cross section PER
     DEUTERON; n(k) is normalised to 1 per deuteron.  Only shapes are
     compared (one free scale per group), so no absolute factor enters.

MEASURED (2026-09-27, this tree): FAIL.  PWIA-window shape chi2/ndf,
187 dof in 57 (Q2, W*, cos) groups of >= 2 points (one single-point group
skipped; 245 of the 604 kept points), stat (+) syst:
  Hulthen (default)  LC 627.85 (p = 5e-49)   instant form 322.54 (p = 3e-9)
  AV18               LC 579.12 (p = 1e-41)   instant form 805.16
  CD-Bonn            LC 278.81 (p = 1.5e-5)  instant form 359.88
  Hulthen LC x GlauberFsiWeight 422.92.
NO wave under EITHER prescription reaches p >= 0.01, so the choice of the
LC prescription for the verdict does not decide it.  Where the default
misses: data / fitted model runs 0.79, 0.95, 1.27, 1.72, 1.41 over the five
p_s bins -- the Hulthen tail falls too fast between 0.30 and 0.46 GeV/c
(AV18 the other way at 0.53: 0.64).  Transverse window: data/PWIA rises
0.76 -> 2.96 from p_s = 0.30 to 0.53; the tree's Glauber weight points the
same way in 4 of 5 bins (not at 0.30) but overshoots 2-6x (1.68 -> 18.9).
CAVEATS on these numbers: the model is evaluated AT the bin-average p_s,
not averaged over the bin (the last bin is wide and n(k) is steep); the
x* off-shell drift above is uncorrected; syst is treated as uncorrelated.

WHAT THIS DOES NOT VALIDATE: the tagged spectrum at EIC kinematics (the
frame caveat of t2_bonus_spectator_shape.py applies: Deeps tags backward in
the target rest frame), p_s < 0.28 GeV/c (BONuS's window, still blocked),
or the absolute rate.

Run:  source env.sh && python3 validation/benchmarks/t2_deeps_spectator_tail.py

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT row -- pass, recorded fail and blocked alike, since the
row carries the verdict; 2 when the harness itself is broken (a missing or
altered vendored file, a missing dependency, any exception).  A harness is
a measurement, not a CI gate: the pytests are the gate.
"""

import collections
import hashlib
import json
import math
import os
import sys
import traceback

HERE = os.path.dirname(os.path.abspath(__file__))
DEEPS = os.path.join(HERE, "data", "clas_deeps_klimenko2006_f2P.json")
#: sha256 of every byte of the vendored file -- the digest
#: t2_bonus_spectator_shape.DEEPS_SHA256 records (2026-09-26, f0a8f1e)
DEEPS_SHA256 = "1d98fd11e1b11b1e070e5504f4919b98b87dbabfe29cdca347522ea4f8737256"

NAME = "t2_deeps_spectator_tail"
#: the paper's two Q2 bins (Klimenko 2006 Sec. IV averages), GeV^2
PAPER_Q2 = (1.8, 2.8)
#: database p_s label (MeV) -> the paper's bin AVERAGE (GeV/c), Sec. IV
PS_PAPER_GEV = {300.0: 0.30, 340.0: 0.34, 390.0: 0.39, 460.0: 0.46,
                560.0: 0.53}
#: PWIA window: cos theta_pq <= this (bin centres); transverse: |cos| < this
PWIA_COS_MAX = -0.3
TRANSVERSE_ABS_COS = 0.3
P_MIN = 0.01            # tolerance: p(chi2, ndf) >= 0.01, declared up front

WAVES = ("hulthen", "av18", "cdbonn")
DEFAULT_WAVE = "hulthen"
WAVE_LABELS = {
    "hulthen": "Hulthen pair, P_D = 0.045 (deuteron_channel() default)",
    "av18": "AV18 (deuteron_channel(source=VmcAV18), fdeut.av18)",
    "cdbonn": "CD-Bonn (cdbonn_wave(), Machleidt 2001 App. D)",
}


# ------------------------------------------------------------------ data

def read_vendored(path=DEEPS, sha256=DEEPS_SHA256):
    """The bytes of the vendored file -- RuntimeError unless their sha256 is
    the recorded one, so a moved number is refused, never silently read."""
    with open(path, "rb") as f:
        raw = f.read()
    got = hashlib.sha256(raw).hexdigest()
    if got != sha256:
        raise RuntimeError("%s: vendored file sha256 %s != recorded %s"
                           % (path, got, sha256))
    return raw


def load_deeps(path=DEEPS):
    return json.loads(read_vendored(path).decode("utf-8"))


def usable_blocks(d):
    """The blocks at the paper's Q2 bins, keyed (q2, W*, p_s label); the key
    is asserted unique (the unverified blocks' labels are not)."""
    out = {}
    for b in d["blocks"]:
        if b["q2"] not in PAPER_Q2:
            continue
        key = (b["q2"], b["w_star"], b["ps_MeV"])
        if key in out:
            raise RuntimeError("duplicate paper-bin block %r" % (key,))
        if b["ps_MeV"] not in PS_PAPER_GEV:
            raise RuntimeError("unknown p_s label %r" % (b["ps_MeV"],))
        out[key] = b
    return out


def census(d):
    blocks = usable_blocks(d)
    rows = [r for b in blocks.values() for r in b["rows"]]
    missing = sum(1 for r in rows if r[1] == 0.0 and r[2] == 0.0)
    return dict(
        blocks_total=len(d["blocks"]),
        blocks_usable=len(blocks),
        blocks_not_read=len(d["blocks"]) - len(blocks),
        q2=sorted({k[0] for k in blocks}),
        w_star=sorted({k[1] for k in blocks}),
        ps_labels_MeV=sorted({k[2] for k in blocks}),
        rows_usable_blocks=len(rows),
        rows_missing_dropped=missing,
        rows_negative_stat0_dropped=sum(1 for r in rows
                                        if r[1] < 0.0 and r[2] == 0.0),
        rows_negative_kept=sum(1 for r in rows if r[1] < 0.0 and r[2] > 0.0),
        rows_kept=sum(1 for r in rows if r[2] > 0.0),
    )


def points(d):
    """Every measured row of the usable blocks as a dict: a row with stat = 0
    (value = 0 -- MISSING -- or one of the two negative values) is dropped."""
    out = []
    for (q2, w, ps_lab), b in sorted(usable_blocks(d).items()):
        for c, val, stat, syst in b["rows"]:
            if stat == 0.0:
                continue            # MISSING (or unweighable), not a zero
            out.append(dict(q2=q2, w=w, ps=PS_PAPER_GEV[ps_lab], c=c,
                            val=val, stat=stat, syst=syst,
                            err=math.hypot(stat, syst)))
    return out


# ------------------------------------------------------------------ tree

class Tree:
    """The tree side: n(k) of each wave, the LC map, the Glauber weight."""

    def __init__(self):
        import lipolgen as lg                   # env.sh puts it on the path
        self.lg = lg
        self.m = lg.M_NUCLEON
        self.m_p = lg.PROTON_MASS
        self.m_d = lg.nuclear_mass(1, 2)
        self.models = {
            "hulthen": lg.TaggedModel(lg.deuteron_channel()),
            "av18": lg.TaggedModel(lg.deuteron_channel(
                source=lg.ClusterWaveSource.VmcAV18)),
        }
        self.cdb = lg.cdbonn_wave()
        self._fsi = None

    def n_of_k(self, wave, k):
        """Unpolarized momentum density, arbitrary normalisation for CD-Bonn
        (shapes only), per deuteron for the TaggedModel waves."""
        if wave == "cdbonn":
            p = k / self.lg.HBARC_GEV_FM
            return self.cdb.psi_s(p) ** 2 + self.cdb.psi_d(p) ** 2
        t = self.models[wave]
        return sum(t.n_of_kc(mi, k, 0.0) for mi in (1.0, 0.0, -1.0)) / 3.0

    def isotropy(self, wave, k):
        """max/min over cos theta_k of the M-averaged density (1 = isotropic)."""
        t = self.models[wave]
        v = [sum(t.n_of_kc(mi, k, c) for mi in (1.0, 0.0, -1.0))
             for c in (-0.9, -0.3, 0.0, 0.4, 0.9)]
        return max(v) / min(v)

    def lc(self, ps, c):
        """(alpha_s, p_T, k_LC, flux E_k / (2 - alpha_s)) for a spectator
        proton of lab momentum ps at cos theta_pq = c."""
        e_s = math.sqrt(self.m_p ** 2 + ps ** 2)
        alpha = (e_s - ps * c) / (0.5 * self.m_d)
        pt2 = ps * ps * max(0.0, 1.0 - c * c)
        k2 = (self.m ** 2 + pt2) / (alpha * (2.0 - alpha)) - self.m ** 2
        k = math.sqrt(max(k2, 0.0))
        return alpha, math.sqrt(pt2), k, math.sqrt(self.m ** 2 + k * k) / (2.0 - alpha)

    def p_model(self, wave, ps, c, form="lc"):
        if form == "instant":
            return self.n_of_k(wave, ps)
        _, _, k, flux = self.lc(ps, c)
        return self.n_of_k(wave, k) * flux

    @property
    def fsi(self):
        if self._fsi is None:
            self._fsi = self.lg.GlauberFsiWeight(self.models[DEFAULT_WAVE])
        return self._fsi

    def fsi_weight(self, ps, c, w, q2):
        lg = self.lg
        x = q2 / (w * w - self.m ** 2 + q2)
        return self.fsi.weight(lg.FsiKinematics(
            k=ps, cos_theta_k=c, phi_k=0.0, w=w, q2=q2, x=x, spectator_z=1,
            spectator_a=1, channel=lg.Channel.TaggedDeuteronP))


# ------------------------------------------------------------------ stats

def chi2_sf(chi2, ndf):
    """Upper-tail probability of chi2 with ndf dof: scipy if present, else
    the T1 harnesses' scipy-free `chi2_sf_fallback` (CI has no scipy)."""
    if ndf <= 0:
        return float("nan")
    try:
        from scipy.stats import chi2 as _c
    except ImportError:
        import importlib.util
        spec = importlib.util.spec_from_file_location(
            "_t1_nmc_f2d", os.path.join(HERE, "t1_nmc_f2d.py"))
        t1 = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(t1)
        return t1.chi2_sf_fallback(chi2, ndf)
    return float(_c.sf(chi2, ndf))


def scaled_chi2(d, m, e):
    """(chi2, ndf, scale) of data d against model m up to ONE free scale:
    s = sum d m / e^2 / sum m^2 / e^2, ndf = n - 1.  Invariant under d -> a d
    (the scale absorbs a) and under m -> b m: a common factor cancels."""
    swm = sum(mi * mi / ei ** 2 for mi, ei in zip(m, e))
    s = sum(di * mi / ei ** 2 for di, mi, ei in zip(d, m, e)) / swm
    chi2 = sum(((di - s * mi) / ei) ** 2 for di, mi, ei in zip(d, m, e))
    return chi2, len(d) - 1, s


def shape_chi2(pts, model, err_key="err"):
    """The PWIA-window shape chi2: one free scale per (Q2, W*, cos) group.
    `model(p)` is the (unnormalised) prediction at point p.  Returns the
    combined, per-Q2 and per-(Q2, W*) sums and the group count."""
    groups = collections.defaultdict(list)
    for p in pts:
        if p["c"] <= PWIA_COS_MAX:
            groups[(p["q2"], p["w"], p["c"])].append(p)
    by_ps = collections.defaultdict(lambda: [0.0, 0.0])
    per_qw = collections.defaultdict(lambda: [0.0, 0])
    per_q = collections.defaultdict(lambda: [0.0, 0])
    tot = [0.0, 0]
    used = skipped = 0
    for (q2, w, _), g in sorted(groups.items()):
        if len(g) < 2:
            skipped += 1
            continue
        used += 1
        mod = [model(p) for p in g]
        c2, ndf, s = scaled_chi2([p["val"] for p in g], mod,
                                 [p[err_key] for p in g])
        for p, mi in zip(g, mod):
            wt = (s * mi / p[err_key]) ** 2     # 1 / variance of the ratio
            by_ps[p["ps"]][0] += wt * p["val"] / (s * mi)
            by_ps[p["ps"]][1] += wt
        for acc in (per_qw[(q2, w)], per_q[q2], tot):
            acc[0] += c2
            acc[1] += ndf
    return dict(
        chi2=tot[0], ndf=tot[1], p=chi2_sf(tot[0], tot[1]),
        groups_used=used, groups_skipped_single_point=skipped,
        # where it misses: the weighted mean of data / (fitted model) per
        # p_s bin, over every group -- 1 everywhere for a perfect shape
        data_over_model_by_ps=[by_ps[k][0] / by_ps[k][1]
                               for k in sorted(by_ps)],
        per_q2={"Q2=%g" % q: dict(chi2=v[0], ndf=v[1], p=chi2_sf(*v))
                for q, v in sorted(per_q.items())},
        per_q2_w={"Q2=%g,W*=%g" % k: dict(chi2=v[0], ndf=v[1], p=chi2_sf(*v))
                  for k, v in sorted(per_qw.items())},
    )


def transverse(pts, tree, wave=DEFAULT_WAVE):
    """Data/PWIA in |cos| < 0.3 against the tree's Glauber direction.
    Recorded only."""
    by_qw = collections.defaultdict(list)
    for p in pts:
        by_qw[(p["q2"], p["w"])].append(p)
    rows = []
    for (q2, w), g in sorted(by_qw.items()):
        bw = [p for p in g if p["c"] <= PWIA_COS_MAX]
        tr = [p for p in g if abs(p["c"]) < TRANSVERSE_ABS_COS]
        if len(bw) < 2 or not tr:
            continue
        m_bw = [tree.p_model(wave, p["ps"], p["c"]) for p in bw]
        f_bw = [tree.fsi_weight(p["ps"], p["c"], w, q2) for p in bw]
        e_bw = [p["err"] for p in bw]
        d_bw = [p["val"] for p in bw]
        _, _, s = scaled_chi2(d_bw, m_bw, e_bw)
        _, _, s_f = scaled_chi2(d_bw, [a * b for a, b in zip(m_bw, f_bw)], e_bw)
        for p in tr:
            m = s * tree.p_model(wave, p["ps"], p["c"])
            rows.append(dict(ps=p["ps"], ratio=p["val"] / m, err=p["err"] / m,
                             tree=(s_f / s) * tree.fsi_weight(p["ps"], p["c"],
                                                              w, q2)))
    per_ps = {}
    for ps in sorted(set(PS_PAPER_GEV.values())):
        sel = [r for r in rows if r["ps"] == ps]
        if not sel:
            continue
        wts = [1.0 / r["err"] ** 2 for r in sel]
        dm = sum(r["ratio"] * wi for r, wi in zip(sel, wts)) / sum(wts)
        tm = sum(r["tree"] * wi for r, wi in zip(sel, wts)) / sum(wts)
        per_ps["ps=%.2f" % ps] = dict(
            n=len(sel), data_over_pwia=dm, data_err=1.0 / math.sqrt(sum(wts)),
            tree_glauber=tm, same_side_of_1=(dm > 1.0) == (tm > 1.0))
    return dict(n_points=len(rows), per_ps=per_ps,
                same_side_all=all(v["same_side_of_1"] for v in per_ps.values()))


# ------------------------------------------------------------------ the row

def measure(path=DEEPS):
    d = load_deeps(path)
    cen = census(d)
    pts = points(d)
    tree = Tree()
    out = dict(census=cen, n_points=len(pts),
               n_points_pwia=sum(1 for p in pts if p["c"] <= PWIA_COS_MAX),
               isotropy_hulthen=tree.isotropy("hulthen", 0.4),
               isotropy_av18=tree.isotropy("av18", 0.4),
               k_lc_range_pwia=[min(tree.lc(p["ps"], p["c"])[2] for p in pts
                                    if p["c"] <= PWIA_COS_MAX),
                                max(tree.lc(p["ps"], p["c"])[2] for p in pts
                                    if p["c"] <= PWIA_COS_MAX)])
    waves = {}
    for wv in WAVES:
        lc = shape_chi2(pts, lambda p, wv=wv: tree.p_model(wv, p["ps"], p["c"]))
        lc_stat = shape_chi2(pts, lambda p, wv=wv: tree.p_model(
            wv, p["ps"], p["c"]), err_key="stat")
        inst = shape_chi2(pts, lambda p, wv=wv: tree.p_model(
            wv, p["ps"], p["c"], form="instant"))
        waves[wv] = dict(lc=lc, lc_stat_only=lc_stat, instant=inst)
    fsi = shape_chi2(pts, lambda p: tree.p_model(DEFAULT_WAVE, p["ps"], p["c"])
                     * tree.fsi_weight(p["ps"], p["c"], p["w"], p["q2"]))
    out.update(waves=waves, lc_times_glauber=fsi,
               transverse=transverse(pts, tree),
               fsi_sigma_xn_mb=tree.fsi.sigma_eff_mb(2.0),
               fsi_clipped_grid_fraction=tree.fsi.clipped_grid_fraction())
    return out


def _sr(r, in_verdict, status=None, **extra):
    if status is None:
        status = "pass" if r["p"] >= P_MIN else "fail"
    return dict(status=status, in_verdict=in_verdict, chi2=r["chi2"],
                ndf=r["ndf"], p=r["p"], groups=r["groups_used"],
                data_over_model_by_ps=r["data_over_model_by_ps"],
                per_q2=r["per_q2"], **extra)


def run(path=DEEPS, verbose=True):
    m = measure(path)
    head = m["waves"][DEFAULT_WAVE]["lc"]
    status = "pass" if head["p"] >= P_MIN else "fail"
    subrows = {}
    for wv in WAVES:
        w = m["waves"][wv]
        subrows["%s/lc" % wv] = _sr(w["lc"], wv == DEFAULT_WAVE,
                                    wave=WAVE_LABELS[wv],
                                    chi2_stat_only=w["lc_stat_only"]["chi2"])
        subrows["%s/instant-form" % wv] = _sr(w["instant"], False,
                                              wave=WAVE_LABELS[wv])
    subrows["hulthen/lc x GlauberFsiWeight"] = _sr(m["lc_times_glauber"], False)
    tr = m["transverse"]
    subrows["transverse window (recorded)"] = dict(
        status="recorded", in_verdict=False, n_points=tr["n_points"],
        same_side_of_1_all=tr["same_side_all"],
        **{k: "data/PWIA %.3f +- %.3f, tree %.3f" % (
            v["data_over_pwia"], v["data_err"], v["tree_glauber"])
           for k, v in tr["per_ps"].items()})
    report = dict(
        name=NAME,
        generator="shape chi2/ndf = %.2f/%d (p = %.3g), Hulthen default, LC "
                  "n(k) E_k/(2-alpha), cos_pq <= -0.3; AV18 %.2f/%d, CD-Bonn "
                  "%.2f/%d" % (head["chi2"], head["ndf"], head["p"],
                               m["waves"]["av18"]["lc"]["chi2"],
                               m["waves"]["av18"]["lc"]["ndf"],
                               m["waves"]["cdbonn"]["lc"]["chi2"],
                               m["waves"]["cdbonn"]["lc"]["ndf"]),
        reference="CLAS Deeps F2N x P(p_s), Klimenko PRC 73 (2006) 035212 "
                  "(60 blocks at Q2 = 1.8/2.8, %d PWIA-window points)"
                  % m["n_points_pwia"],
        tolerance="p(chi2, ndf) >= %.2f, default wave, one free scale per "
                  "(Q2, W*, cos) group" % P_MIN,
        status=status,
        config=dict(wave=WAVE_LABELS[DEFAULT_WAVE],
                    prescription="LC (Frankfurt-Strikman), coded in harness",
                    pwia_cos_max=PWIA_COS_MAX,
                    transverse_abs_cos=TRANSVERSE_ABS_COS,
                    errors="stat (+) syst in quadrature",
                    fsi_sigma_xn_mb=m["fsi_sigma_xn_mb"]),
        subrows=subrows,
        measured=m,
    )
    if verbose:
        c = m["census"]
        print("Deeps file: %d blocks, %d read (Q2 = %s), %d not read; %d rows, "
              "%d missing (value = stat = 0) and %d negative with stat = 0 "
              "dropped, %d kept"
              % (c["blocks_total"], c["blocks_usable"], c["q2"],
                 c["blocks_not_read"], c["rows_usable_blocks"],
                 c["rows_missing_dropped"], c["rows_negative_stat0_dropped"],
                 c["rows_kept"]))
        print("LC k over the PWIA-window points: %.3f .. %.3f GeV; M-averaged "
              "n(k, c) max/min over c at k = 0.4: Hulthen %.12f, AV18 %.12f"
              % (tuple(m["k_lc_range_pwia"]) + (m["isotropy_hulthen"],
                                                m["isotropy_av18"])))
        print("%-32s %14s %10s %18s %14s" % ("PWIA-window shape", "chi2/ndf",
                                              "p", "stat-only chi2", "verdict"))
        for name, s in subrows.items():
            if "chi2" not in s:
                continue
            print("%-32s %9.2f/%-4d %10.3g %18s %14s" % (
                name, s["chi2"], s["ndf"], s["p"],
                "%.2f" % s["chi2_stat_only"] if "chi2_stat_only" in s else "",
                s["status"].upper() if s["in_verdict"] else "(recorded)"))
        print("data / fitted model per p_s bin %s GeV/c:"
              % [round(v, 2) for v in sorted(set(PS_PAPER_GEV.values()))])
        for name, s in subrows.items():
            if "data_over_model_by_ps" in s:
                print("  %-30s %s" % (name, " ".join(
                    "%.3f" % v for v in s["data_over_model_by_ps"])))
        print("default wave, per (Q2, W*):")
        for k, v in head["per_q2_w"].items():
            print("  %-16s %7.2f/%-3d p = %.3g" % (k, v["chi2"], v["ndf"], v["p"]))
        print("transverse |cos_pq| < 0.3 (recorded, not asserted), %d points; "
              "GlauberFsiWeight sigma_XN = %.0f mb, clipped grid %.4f:"
              % (tr["n_points"], m["fsi_sigma_xn_mb"],
                 m["fsi_clipped_grid_fraction"]))
        for k, v in tr["per_ps"].items():
            print("  %s: data/PWIA %.3f +- %.3f, tree Glauber %.3f, same side "
                  "of 1: %s" % (k, v["data_over_pwia"], v["data_err"],
                                v["tree_glauber"], v["same_side_of_1"]))
        print("REPORT | %s | %s | %s | %s | %s" % (
            NAME, report["generator"], report["reference"],
            report["tolerance"], status.upper()))
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
