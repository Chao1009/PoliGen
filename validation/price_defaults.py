#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""price_defaults -- the default-input price table of the RANKED observables.

Run 2026-10-26, phase A3 of `docs/open_items/run_2026-10-26/PLAN.md`; asked
for by `docs/benchmarking/08_channel_priority.md` sec. 5.1.  The record that
carries the table and reads it is
`docs/open_items/run_2026-10-26/phase_A_baseline.md`, section "A3 -- the
price table".

A RECIPE, NOT A BENCHMARK HARNESS.  The name does not match `t[1-6]_*.py`, so
`make_report.py` never runs it: there is no reference, no tolerance and no
verdict.  It PRICES.  For each channel it evaluates the channel's OWN ranked
observable under each knob that the knob-provenance table
(`meta["knob_provenance"]`, docs/USAGE.md sec. 7c) reports as READ on that
channel and plan, and quotes it as the ratio to the shipped default on the
same channel, plan, window, N and seed.  A knob the table reports `not-read`
is recorded by its LABEL, and the run is checked bit-identical to the default
one; a knob `PipelineConfig.validate()` refuses is recorded by its REFUSAL
TEXT; a backend whose build tier or LHAPDF set is missing is recorded as
`environment` with the library's own message.  None of the three is ever
written as a number.  NOTHING HERE MOVES A DEFAULT.

THE WINDOW -- every number is in it unless its row says otherwise
  beam      default_configs(<isotope>)[1]: e 10 GeV x ion 99.5 GeV/u
  scenario  the shipped one (x_max = 1.0); the `li6-convolution` rows run at
            x_max = 0.95, which that backend needs (docs/USAGE.md sec. 2a,
            "The top x cell"); the four super-bins sit at x <= 0.2 and their
            cells do not move with it
  optics    yr-high-acceptance (the lg.run default; it decides the tagged
            channels' Roman-Pot route and the A_zz^tag acceptance weights)
  run       lg.run(events=N, seed=SEED), N = 200 000, SEED = 20261026, in
            fixed-count mode; the plan fills at lg.run's own defaults
            P_z = 0.7, P_zz = 0.6, P_e = 0.7 (P_e is read by helicity-flip
            only), pzz_mode = ladder

THE OBSERVABLES, per channel (the rank is 08_channel_priority.md sec. 2)
  inclusive 6Li, rank 1 (claim 1) -- A2(S): the cos 2phi' amplitude per unit
      P_zz, a2_eff / P_zz, on the `transverse-tensor` fill (theta_S = pi/2,
      the plan's one category), in the four sweet-spot SUPER-BINS of the PLB
      letter (`SWEET_SPOTS`).  a2_eff is `InclusiveSampler::
      effective_modulation` (src/core/sampler.cpp) restricted to the run's
      own accepted cells whose centres lie in the super-bin -- the
      restriction the C++ takes as `mask` and the binding does not expose,
      re-done here from `state_tables` and checked against the unmasked
      library call to 1e-12.  This is the sibling's own estimator target
      (PolarizedLithiumSim/evgen/scripts/money_cos2phi.py `measure`, 'truth'
      = a2_eff / pzz; `superbin_mask` selects by cell centre).
  inclusive 6Li, rank 4 (paper 2) -- the run-level A_zz: the canonical
      thirds estimator (`estimators.azz_thirds`, which divides by P_zz) on
      the `tensor-thirds` per-category cross sections in meta.
  inclusive 7Li, rank 5 (paper 3) -- A_par = `estimators.apar_flip` on the
      `helicity-flip` per-category cross sections in meta, P_e = P_z = 0.7.
  tagged 7Li alpha, rank 3 -- A_par as above, and the spectator spectrum's
      mean k from the generated `k` column (standard error std/sqrt(N)).
  tagged 6Li alpha, rank 6 -- A_zz^tag(k): `azz_tensor_curve_weighted` of the
      run's own TaggedModel with `acceptance_weights` at the run's own optics
      (n_phi = 32, the setting of 07_cw_sign_investigation.md sec. 6.4), at
      the grid cells nearest k = 0.20 and 0.30 GeV.  Cross-checked against
      the GENERATED Roman-Pot-accepted events, thirds estimator in
      k +- 0.025 GeV, against the model's acceptance-weighted average over
      the same bin.  With no acceptance weight the same curve is ~1e-4 (the
      family's note prints it): the S-D tensor structure is a cos theta_k
      shape, so A_zz^tag exists only through the acceptance.
  coherent 6Li, rank 2 (claim 2) -- c2 of the recoil weight
      1 + c2 cos 2(phi_t - phi_S) on the `transverse-tensor` fill: the
      GENERATED modulation 2 sum(w cos 2(phi_t - phi_S)) / sum(w) (the
      default azimuth is uniform with the modulation in `weight`,
      src/core/coherent.cpp), its model expectation over the same |t|
      sample, and the model c2 at the tag <|t|> of the letter's 0.20 GeV
      envelope (`CoherentScenario::mean_t_tagged`; no coherent recoil reaches
      the pots at the YR envelope in this sample).
  every channel -- sigma (meta["sigma_pb"]) as the rate line; on the tagged
      channels also the Roman-Pot tag fraction of the generated sample.

WHAT IS EXACT AND WHAT IS SAMPLED.  sigma, the per-category cross sections,
A2(S), A_zz, A_par, A_zz^tag(k) and the model c2 are N-independent functions
of the configuration.  <k>, the tag fractions, the generated c2 and every
"generated" cross-check are sampled at N, and carry their errors.  In
fixed-count mode the per-category event counts are a largest-remainder
ROUNDING of N in proportion to the exact cross sections (src/core/
pipeline.cpp, the `counts` block), so a counting A_zz or A_par on those
columns carries rounding, not statistics: the prices use meta.

Run:   source env.sh && python3 validation/price_defaults.py [--plot]
       [--events N] [--seed S] [--families li6-incl,li7-incl,...]
       [--json PATH] [--out-dir DIR]
`--plot` draws docs/benchmarking/figures/price_defaults.png through
validation/benchmarks/_plots.py (or DIR/price_defaults.png).  Exit 0 when the
table was printed, 2 on any exception (the harness convention, kept although
this is not a harness).
"""

import argparse
import hashlib
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
BENCH = os.path.join(HERE, "benchmarks")

NAME = "price_defaults"
N_DEFAULT = 200000
SEED = 20261026
PZ, PZZ, PE = 0.7, 0.6, 0.7       # lg.run's own plan defaults
X_MAX_CONVOLUTION = 0.95          # docs/USAGE.md sec. 2a, "The top x cell"
N_PHI_ACCEPTANCE = 32             # 07_cw_sign_investigation.md sec. 6.4
K_POINTS = (0.20, 0.30)           # GeV; read at the nearest TaggedModel cell
K_HALF_WIDTH = 0.025              # GeV; the generated cross-check's k bin
TAG_ENVELOPE_PT = 0.20            # GeV; the letter's coherent envelope
EPS_B0_MEASURED = -0.0070         # coherent.hpp, eps_b0: the measured-Q row
F0_BAND = (0.02, 0.08)            # coherent.hpp, f0: the scenario band
RP_ROUTES = (1, 4)                # spectator.hpp `rp_accepted`: RP main + near-beam
#: shown, never priced: c2_gen is sampled (its ratio is MC noise) and the
#: tagged-6Li sigma_tot A_zz is zero by construction (docs/USAGE.md sec. 2b)
NO_RATIO = ("c2_gen", "azz_sigma_tot")

#: The four sweet-spot SUPER-BINS of the letter's claim 1, copied from the
#: sibling repository (read only): PolarizedLithiumSim/paper/figs/fig1.json,
#: key "superbins", read 2026-10-26.  Each is the 40 x 30 analysis-map cell of
#: highest cos 2phi significance in its Q2 band (money_cos2phi.py
#: `pick_sweet_spots_banded`) padded by one bin on every side
#: (`superbin_edges(pad=1)`, clipped at the Q2 = 1 floor); plans/07 sec. 7.1
#: lists the centres as spot 1 (x 0.0282, Q2 1.14), spot 2 (0.0112, 1.14),
#: spot 3 (0.0708, 3.13), spot 4 (0.141, 14.3).
#: (spot, x centre, Q2 centre [GeV^2], x_lo, x_hi, Q2_lo, Q2_hi)
SWEET_SPOTS = (
    (1, 0.02818382931264455, 1.1350556807493803,
     0.019952623149688806, 0.039810717055349734, 1.0, 1.6598493257629445),
    (2, 0.011220184543019636, 1.1350556807493803,
     0.00794328234724282, 0.01584893192461114, 1.0, 1.6598493257629445),
    (3, 0.07079457843841383, 3.127191661128126,
     0.05011872336272725, 0.09999999999999999,
     2.1384691999823757, 4.573050519273263),
    (4, 0.14125375446227548, 14.300805449788998,
     0.09999999999999999, 0.19952623149688806,
     9.779327685429287, 20.91279105182546),
)

#: knob -> its row name in meta["knob_provenance"]
KP_NAME = {"unpol_sf": "unpol_sf", "pol_sf": "pol_sf", "r_source": "r_source",
           "cluster_wave": "cluster_wave", "triton_sf": "triton_sf",
           "eps_b0": "coherent_eps_b0", "f0": "coherent_f0",
           "b1_unpol": "b1_unpol", "b1_band_scale": "b1_band_scale"}

FAMILIES = ("li6-incl", "li6-incl-r", "li7-incl", "li7-tag", "li6-tag",
            "li6-coh")

_ENV_WORDS = ("LHAPDF", "PYTHIA", "built without", "not installed",
              "Info file not found")


def _import_lipolgen():
    try:
        import lipolgen
    except ImportError:
        sys.path.insert(0, os.path.join(ROOT, "build", "python"))
        import lipolgen
    return lipolgen


lg = _import_lipolgen()


# ------------------------------------------------------------------ helpers

def _attempt(kw, events, seed):
    """`lg.run(events, seed, **kw)` -> (out, None), or (None, message) when
    the library refuses the configuration or misses a tier (RuntimeError)."""
    try:
        return lg.run(events=events, seed=seed, **kw), None
    except RuntimeError as e:
        return None, str(e)


def _is_environment(msg):
    """A missing build tier or LHAPDF set, as the library words it -- never a
    `PipelineConfig: ...` refusal, which is a statement about the knob."""
    return (not msg.startswith("PipelineConfig")
            and any(w in msg for w in _ENV_WORDS))


def _kp(out, knob):
    row = out["meta"]["knob_provenance"][KP_NAME.get(knob, knob)]
    return dict(status=row["status"], value=row.get("value"),
                label=row.get("label") or "", reason=row.get("reason") or "")


def _fingerprint(out):
    """sha256 over the generated columns every channel writes and over the
    meta cross sections: equal fingerprints = the same file, bit for bit."""
    h = hashlib.sha256()
    for c in ("x", "q2", "y", "phi", "weight", "k", "t", "phi_spec", "route",
              "category_index"):
        if c in out:
            h.update(np.ascontiguousarray(out[c]).tobytes())
    m = out["meta"]
    h.update(np.asarray([m["sigma_pb"]] + list(m["sigma_per_category_pb"]),
                        dtype=float).tobytes())
    return h.hexdigest()


def _rp_fraction(out):
    acc = np.isin(np.asarray(out["route"]), RP_ROUTES)
    n = acc.size
    f = float(acc.mean())
    return f, math.sqrt(max(f * (1.0 - f), 0.0) / n), int(acc.sum())


def _row(setting, knob=None, value=None, provenance="read", label="",
         values=None, errors=None, identical=None):
    return dict(setting=setting, knob=knob, value=value, provenance=provenance,
                label=label, values=dict(values or {}),
                errors=dict(errors or {}), ratios={}, identical=identical)


def _probe(base_kw, knob, value, events, seed, default_row, extract,
           expect):
    """A knob the table does not price: `expect` is "not-read" (run it, take
    the label, check the file is bit-identical to the default) or "refused"
    (validate() must throw; keep its text)."""
    kw = dict(base_kw)
    if knob in ("eps_b0", "f0"):
        kw["coherent"] = {knob: value}
    else:
        kw[knob] = value
    setting = "%s = %s" % (knob, value)
    out, err = _attempt(kw, events, seed)
    if err is not None:
        if _is_environment(err):
            return _row(setting, knob, value, "environment", err)
        if expect != "refused":
            raise RuntimeError("%s: expected %s, refused: %s"
                               % (setting, expect, err))
        return _row(setting, knob, value, "refused", err)
    kp = _kp(out, knob)
    if expect == "refused":
        raise RuntimeError("%s: expected a refusal, ran with status %s"
                           % (setting, kp["status"]))
    if kp["status"] != "not-read":
        raise RuntimeError("%s: expected not-read, provenance says %s"
                           % (setting, kp["status"]))
    vals, _errs = extract(out)
    same = (_fingerprint(out) == default_row["fingerprint"]
            and all(default_row["values"].get(k) == v for k, v in vals.items()))
    r = _row(setting, knob, value, "not-read",
             kp["label"] or kp["reason"].split(":")[0], identical=same)
    del out
    return r


def _ratios(fam):
    """value / default value, column by column, for every row that ran."""
    base = fam["rows"][0]["values"]
    for r in fam["rows"][1:]:
        if r["provenance"] not in ("read",):
            continue
        for k, v in r["values"].items():
            b = base.get(k)
            if k in NO_RATIO:
                continue
            if isinstance(v, (int, float)) and isinstance(b, (int, float)) \
                    and b != 0 and math.isfinite(v) and math.isfinite(b):
                r["ratios"][k] = v / b


# ------------------------------------------------------- the observables

def _superbin_amplitudes(out):
    """A2(S) = a2_eff/P_zz in each super-bin, from the run's own sampler
    tables and its plan's one category; plus the generated-column moment
    estimator on the events whose cell lies in the super-bin."""
    p = out["pipeline"]
    s = p.dis_sampler
    cats = list(p.plan.categories)
    assert len(cats) == 1, "transverse-tensor has one category"
    cat = cats[0]
    pzz = cat.moments().tensor
    ms = (1.0, 0.0, -1.0)
    tabs = [(pm, s.state_tables(cat, m)) for pm, m in zip(cat.populations, ms)
            if pm > 0.0]
    xs = np.asarray(s.cell_xsec_pb)
    xc, qc = np.asarray(s.x_cells), np.asarray(s.q2_cells)
    one = np.ones_like(xs, dtype=bool)

    def eff(mask):
        den = sum(pm * np.sum(xs[mask] * (1.0 + t["w_avg"][mask])) for pm, t in tabs)
        num = sum(pm * np.sum(xs[mask] * t["a2"][mask]) for pm, t in tabs)
        return den, num / den

    # the restriction must reproduce the library's unmasked call
    lib = s.effective_modulation(cat)
    den_all, a2_all = eff(one)
    if abs(den_all / lib.sigma_pb - 1.0) > 1e-12 or abs(a2_all - lib.a2) > 1e-12 * abs(lib.a2):
        raise RuntimeError("super-bin restriction does not reproduce "
                           "effective_modulation (%r vs %r)" % ((den_all, a2_all), (lib.sigma_pb, lib.a2)))
    mix_a2 = sum(pm * t["a2"] for pm, t in tabs)
    cell = np.asarray(out["cell"])
    php = np.asarray(out["phi"]) - np.asarray(out["phi_s"])
    vals, errs = {}, {}
    for (i, _xc, _qc, xlo, xhi, qlo, qhi) in SWEET_SPOTS:
        mask = (xc >= xlo) & (xc < xhi) & (qc >= qlo) & (qc < qhi)
        _den, a2 = eff(mask)
        vals["A2_S%d" % i] = a2 / pzz
        vals["cells_S%d" % i] = int(mask.sum())
        sel = mask[cell]
        n = int(sel.sum())
        vals["n_S%d" % i] = n
        if n > 1:
            c = np.cos(2.0 * php[sel])
            vals["gen_A2_S%d" % i] = float(lg.estimators.cos2phi_moment(list(php[sel]), pzz))
            errs["gen_A2_S%d" % i] = 2.0 * float(np.std(c)) / math.sqrt(n) / pzz
    vals["sigma_tt"] = out["meta"]["sigma_pb"]
    return vals, errs, (xc, qc, mix_a2)


def _thirds(out):
    m = out["meta"]
    names = list(out["category_names"])
    assert names == ["azz+", "azz-", "azz0"], names
    s = m["sigma_per_category_pb"]
    return dict(sigma=m["sigma_pb"],
                azz=lg.estimators.azz_thirds(s[0], s[1], s[2], PZZ),
                sig_cat=list(s))


def _apar(out):
    m = out["meta"]
    names = list(out["category_names"])
    assert names == ["apar+", "apar-"], names
    s = m["sigma_per_category_pb"]
    return dict(sigma=m["sigma_pb"],
                apar=lg.estimators.apar_flip(s[0], s[1], PE, PZ),
                sig_cat=list(s))


def _tagged_li7(out):
    v = _apar(out)
    k = np.asarray(out["k"])
    f, ef, n_acc = _rp_fraction(out)
    v.update(k_mean=float(k.mean()), rp_frac=f, rp_n=n_acc)
    return v, dict(k_mean=float(k.std()) / math.sqrt(k.size), rp_frac=ef)


def _azz_err(n_p, n_m, n_0):
    s = n_p + n_m + n_0
    if s <= 0 or n_0 <= 0 or n_p + n_m <= 0:
        return float("nan")
    return 3.0 * math.sqrt(n_0 * (n_p + n_m) / s ** 3) / PZZ


def _tagged_li6(out):
    p = out["pipeline"]
    tm = p.tagged_model
    k = np.asarray(tm.k)
    p_u = p.beam_config.ion_momentum_per_nucleon
    eps = np.asarray(lg.acceptance_weights(tm, p_u, p.optics, p.pot_config,
                                           N_PHI_ACCEPTANCE))
    curve = np.asarray(lg.azz_tensor_curve_weighted(tm, eps))
    # the same curve with no acceptance at all (every cos theta_k weight 1)
    curve_4pi = np.asarray(lg.azz_tensor_curve_weighted(tm, np.ones_like(eps)))
    n_p = np.asarray(tm.n_of_kc_table(1.0))
    n_0 = np.asarray(tm.n_of_kc_table(0.0))
    n_m = np.asarray(tm.n_of_kc_table(-1.0))
    # n_M is a density with  int n k^2 dk dOmega = 1  (tagged.hpp): a k-bin
    # average weights each cell by k^2 (uniform dk and dc cancel)
    num_k = ((n_p + n_m - 2.0 * n_0) * eps).sum(axis=1) * k ** 2
    den_k = ((n_p + n_m + n_0) * eps).sum(axis=1) * k ** 2
    names = list(out["category_names"])
    assert names == ["azz+", "azz-", "azz0"], names
    acc = np.isin(np.asarray(out["route"]), RP_ROUTES)
    kg = np.asarray(out["k"])
    ci = np.asarray(out["category_index"])
    t3 = _thirds(out)
    # A_zz on sigma_tot is zero BY CONSTRUCTION here (the struck cluster has no
    # b1 at --inclusive-b1 off): kept as evidence, never as a price
    v = dict(sigma=t3["sigma"], azz_sigma_tot=t3["azz"], sig_cat=t3["sig_cat"])
    e = {}
    for kk in K_POINTS:
        tag = "%03d" % round(kk * 100)
        i = int(np.argmin(np.abs(k - kk)))
        v["k_cell_" + tag] = float(k[i])
        v["azz_tag_" + tag] = float(curve[i])
        v["azz_4pi_" + tag] = float(curve_4pi[i])
        lo, hi = kk - K_HALF_WIDTH, kk + K_HALF_WIDTH
        mi = (k >= lo) & (k < hi)
        v["azz_tag_bin_" + tag] = float(num_k[mi].sum() / den_k[mi].sum())
        sel = acc & (kg >= lo) & (kg < hi)
        n = [int(np.sum(sel & (ci == j))) for j in range(3)]
        v["gen_n_" + tag] = n
        if sum(n) > 0:
            v["gen_azz_" + tag] = float(lg.estimators.azz_thirds(
                float(n[0]), float(n[1]), float(n[2]), PZZ))
            e["gen_azz_" + tag] = _azz_err(*n)
    f, ef, n_acc = _rp_fraction(out)
    v.update(rp_frac=f, rp_n=n_acc, k_mean=float(kg.mean()),
             vec_dil=float(tm.vector_dilution()), tens_dil=float(tm.tensor_dilution()))
    e.update(rp_frac=ef, k_mean=float(kg.std()) / math.sqrt(kg.size))
    return v, e


def _coherent(out):
    p = out["pipeline"]
    sc = p.config.coherent
    names = list(out["category_names"])
    assert len(names) == 1, names
    pzz = p.plan.categories[0].moments().tensor
    t = np.asarray(out["t"])
    w = np.asarray(out["weight"])
    c = np.cos(2.0 * (np.asarray(out["phi_spec"]) - np.asarray(out["phi_s"])))
    a, b = w * c, w
    n = w.size
    r = a.mean() / b.mean()
    var = (np.var(a) - 2.0 * r * np.mean((a - a.mean()) * (b - b.mean()))
           + r * r * np.var(b)) / (n * b.mean() ** 2)
    model = np.array([sc.cos2phi_coefficient(float(ti), pzz) for ti in t])
    # c2 is LINEAR in |t| (coherent.hpp `cos2phi_coefficient`), so its
    # expectation over the generated spectrum is c2 at the spectrum's mean:
    # B exp(-B|t|) renormalised on [0, t_max] (`sample_t`) has
    # <|t|> = 1/B - t_max exp(-B t_max) / (1 - exp(-B t_max)) -- N-independent
    b, tmax = float(sc.slope_b), float(p.config.coherent_t_max)
    t_exp = 1.0 / b - tmax * math.exp(-b * tmax) / (1.0 - math.exp(-b * tmax))
    t_tag = sc.mean_t_tagged(TAG_ENVELOPE_PT)
    _f, _ef, n_acc = _rp_fraction(out)
    v = dict(sigma=out["meta"]["sigma_pb"], c2_gen=float(2.0 * r),
             c2_sample=float(model.mean()), c2_exp=float(sc.cos2phi_coefficient(t_exp, pzz)),
             c2_tag=float(sc.cos2phi_coefficient(t_tag, pzz)), t_exp=float(t_exp),
             t_tag=float(t_tag), t_mean=float(t.mean()), rp_n=n_acc,
             eps_b0=float(sc.eps_b0), f0=float(sc.f0), pzz=float(pzz),
             slope_b=b, t_max=tmax)
    return v, dict(c2_gen=float(2.0 * math.sqrt(max(var, 0.0))),
                   t_mean=float(t.std()) / math.sqrt(n))


# ------------------------------------------------------------- families

def _family(name, title, window, columns):
    return dict(name=name, title=title, window=window, columns=columns,
                rows=[], notes=[])


def _price(fam, setting, knob, value, kw, extract, events, seed, log):
    out, err = _attempt(kw, events, seed)
    if err is not None:
        if _is_environment(err):
            fam["rows"].append(_row(setting, knob, value, "environment", err))
            return None
        raise RuntimeError("%s: refused where it should run: %s" % (setting, err))
    if knob is not None:
        kp = _kp(out, knob)
        if kp["status"] != "read":
            raise RuntimeError("%s: priced but provenance says %s (%s)"
                               % (setting, kp["status"], kp["label"]))
    vals, errs = extract(out)
    r = _row(setting, knob, value, "read" if knob else "default",
             values=vals, errors=errs)
    r["fingerprint"] = _fingerprint(out)
    fam["rows"].append(r)
    log("  %-48s %s" % (fam["name"] + ": " + setting, "ok"))
    return out


def family_li6_inclusive(events, seed, log):
    fam = _family(
        "li6-incl",
        "Inclusive 6Li -- A2(S), the cos 2phi' amplitude of claim 1 (rank 1), "
        "and the run-level A_zz of paper 2 (rank 4)",
        "6Li, channel inclusive; A2 on plan transverse-tensor (P_zz = 0.6, "
        "theta_S = pi/2), A_zz and sigma on plan tensor-thirds (P_z = 0.7, "
        "P_zz = 0.6); x_max 1.0; N = %d, seed %d" % (events, seed),
        [("sigma", "sigma [pb] (tensor-thirds)", "pb"),
         ("azz", "A_zz (run level)", "asym"),
         ("A2_S1", "A2(S1) x 0.028, Q2 1.14", "asym"),
         ("A2_S2", "A2(S2) x 0.011, Q2 1.14", "asym"),
         ("A2_S3", "A2(S3) x 0.071, Q2 3.13", "asym"),
         ("A2_S4", "A2(S4) x 0.141, Q2 14.3", "asym")])
    base_tt = dict(isotope="6Li", channel="inclusive", plan="transverse-tensor")
    base_3 = dict(isotope="6Li", channel="inclusive", plan="tensor-thirds")
    cells = None

    def both(kw, setting, knob, value):
        nonlocal cells
        tt, err = _attempt(dict(base_tt, **kw), events, seed)
        if err is not None:
            if _is_environment(err):
                fam["rows"].append(_row(setting, knob, value, "environment", err))
                return
            raise RuntimeError("%s refused: %s" % (setting, err))
        th, err = _attempt(dict(base_3, **kw), events, seed)
        if err is not None:
            raise RuntimeError("%s refused on tensor-thirds: %s" % (setting, err))
        for o in ((tt, th) if knob else ()):
            kp = _kp(o, knob)
            if kp["status"] != "read":
                raise RuntimeError("%s: priced but provenance says %s" % (setting, kp["status"]))
        v, e, (xc, qc, mix) = _superbin_amplitudes(tt)
        t3 = _thirds(th)
        v.update(sigma=t3["sigma"], azz=t3["azz"])
        r = _row(setting, knob, value, "read" if knob else "default", values=v, errors=e)
        r["fingerprint"] = _fingerprint(th)
        r["fingerprint_tt"] = _fingerprint(tt)
        if cells is None:
            cells = (xc, qc, mix)
        else:
            # F2 cancels in a2 cell by cell: measure it on the common grid
            ok = (len(xc) == len(cells[0]) and np.array_equal(xc, cells[0])
                  and np.array_equal(qc, cells[1]))
            nz = cells[2] != 0.0
            r["values"]["cell_max_rel_a2"] = (
                float(np.max(np.abs(mix[nz] / cells[2][nz] - 1.0))) if ok else float("nan"))
        fam["rows"].append(r)
        log("  %-48s ok" % ("li6-incl: " + setting))
        del tt, th

    both({}, "unpol_sf = toy (shipped default)", None, "toy")
    both(dict(unpol_sf="mstw"), "unpol_sf = mstw", "unpol_sf", "mstw")
    both(dict(unpol_sf="ct18nlo"), "unpol_sf = ct18nlo", "unpol_sf", "ct18nlo")
    default = fam["rows"][0]

    def amp_extract(out):
        v, _e, _c = _superbin_amplitudes(out)
        return {k: v[k] for k in v if k.startswith("A2_")}, {}

    def thirds_extract(out):
        t = _thirds(out)
        return dict(sigma=t["sigma"], azz=t["azz"]), {}

    # pol_sf: NOT READ under any tensor fill -- label, and bit-identity on BOTH plans
    for plan, ext, fp in (("transverse-tensor", amp_extract, "fingerprint_tt"),
                          ("tensor-thirds", thirds_extract, "fingerprint")):
        d = dict(default, fingerprint=default[fp])
        r = _probe(dict(isotope="6Li", channel="inclusive", plan=plan),
                   "pol_sf", "nnpdfpol", events, seed, d, ext, "not-read")
        r["setting"] += " (%s)" % plan
        fam["rows"].append(r)
    d3 = dict(default)
    for knob, val, expect in (("r_source", "r1998", "refused"),
                              ("cluster_wave", "vmc", "refused")):
        fam["rows"].append(_probe(base_3, knob, val, events, seed, d3,
                                  thirds_extract, expect))
    dtt = dict(default, fingerprint=default["fingerprint_tt"])
    for knob, val in (("triton_sf", "ciofi-simula"), ("eps_b0", EPS_B0_MEASURED),
                      ("f0", F0_BAND[0])):
        fam["rows"].append(_probe(base_tt, knob, val, events, seed, dtt,
                                  amp_extract, "not-read"))
    _ratios(fam)
    fam["notes"].append(
        "A2 per cell is the same double on every --unpol-sf backend (F2 cancels: "
        "max |a2_cell/a2_cell(toy) - 1| = %s over the %d accepted cells); the "
        "super-bin A2 moves only through the cross-section weights sigma_c (1 + w_avg) "
        "of its cells." % (", ".join("%.2g" % r["values"].get("cell_max_rel_a2", float("nan"))
                                     for r in fam["rows"][1:3] if r["provenance"] == "read"),
                           len(cells[0]) if cells else 0))
    return fam


def family_li6_inclusive_r(events, seed, log):
    """--r-source on the one path that accepts it."""
    fam = _family(
        "li6-incl-r",
        "Inclusive 6Li on the opt-in li6-convolution path -- --r-source "
        "(refused on the shipped miller b1, so priced here and labelled)",
        "6Li, channel inclusive, --b1-model li6-convolution, --b1-unpol toy, "
        "--b1-band-scale 1 unless stated, x_max %.2f (the backend's top-cell "
        "positivity); A2 on transverse-tensor, A_zz and sigma on tensor-thirds; "
        "N = %d, seed %d" % (X_MAX_CONVOLUTION, events, seed),
        [("sigma", "sigma [pb] (tensor-thirds)", "pb"),
         ("azz", "A_zz (run level, band 1)", "asym"),
         ("A2_S1", "A2(S1)", "asym"), ("A2_S2", "A2(S2)", "asym"),
         ("A2_S3", "A2(S3)", "asym"), ("A2_S4", "A2(S4)", "asym")])
    sc = lg.PipelineConfig().scenario
    sc.x_max = X_MAX_CONVOLUTION
    conv = dict(isotope="6Li", channel="inclusive", scenario=sc,
                b1_model="li6-convolution")
    kernels = {}
    for setting, knob, value in (("r_source = unset (the path's default)", None, "unset"),
                                 ("r_source = r1998", "r_source", "r1998")):
        kw = {} if knob is None else {knob: value}
        tt, err = _attempt(dict(conv, plan="transverse-tensor", **kw), events, seed)
        if err is not None:
            if _is_environment(err):
                fam["rows"].append(_row(setting, knob, value, "environment", err))
                continue
            raise RuntimeError("%s: %s" % (setting, err))
        th, err = _attempt(dict(conv, plan="tensor-thirds", **kw), events, seed)
        if err is not None:
            raise RuntimeError("%s (tensor-thirds): %s" % (setting, err))
        if knob:
            for o in (tt, th):
                if _kp(o, knob)["status"] != "read":
                    raise RuntimeError("%s: not read" % setting)
        v, e, _c = _superbin_amplitudes(tt)
        t3 = _thirds(th)
        v.update(sigma=t3["sigma"], azz=t3["azz"])
        r = _row(setting, knob, value, "read" if knob else "default", values=v, errors=e)
        r["fingerprint"] = _fingerprint(th)
        fam["rows"].append(r)
        kernels[value] = th["pipeline"].dis_sampler.kernel
        log("  %-48s ok" % ("li6-incl-r: " + setting))
        del tt, th
    if not fam["rows"] or fam["rows"][0]["provenance"] != "default":
        return fam
    default = fam["rows"][0]

    def thirds_extract(out):
        t = _thirds(out)
        return dict(sigma=t["sigma"], azz=t["azz"]), {}

    r = _probe(dict(conv, plan="tensor-thirds"), "r_source", "sigma-lt",
               events, seed, default, thirds_extract, "not-read")
    fam["rows"].insert(1, r)
    _ratios(fam)
    # the band: the ratio does not depend on the mandatory band scale
    band = {}
    for scale in (2.0, 0.0):
        for value in ("unset", "r1998"):
            kw = dict(conv, plan="tensor-thirds", b1_band_scale=scale)
            if value != "unset":
                kw["r_source"] = value
            out, err = _attempt(kw, events, seed)
            if err is not None:
                band[(scale, value)] = err
                continue
            band[(scale, value)] = _thirds(out)["azz"]
            del out
    b2u, b2r = band.get((2.0, "unset")), band.get((2.0, "r1998"))
    b0u, b0r = band.get((0.0, "unset")), band.get((0.0, "r1998"))
    if all(isinstance(v, float) for v in (b2u, b2r, b0u, b0r)):
        fam["notes"].append(
            "band: at --b1-band-scale 2, A_zz = %.7e (unset) / %.7e (r1998), ratio "
            "%.6f; at band 0, A_zz = %.3e / %.3e (zero on both: no ratio)." % (
                b2u, b2r, b2r / b2u, b0u, b0r))
        fam["band"] = dict(azz_unset_2=b2u, azz_r1998_2=b2r, ratio_2=b2r / b2u,
                           azz_unset_0=b0u, azz_r1998_0=b0r)
    else:
        fam["notes"].append("band check not complete: %r" % band)
    # the quotable configuration of the A = 2 gate: --b1-unpol mstw
    mst = {}
    for value in ("unset", "r1998"):
        kw = dict(conv, plan="tensor-thirds", b1_unpol="mstw")
        if value != "unset":
            kw["r_source"] = value
        out, err = _attempt(kw, events, seed)
        mst[value] = err if err is not None else _thirds(out)
        out = None
    if all(isinstance(v, dict) for v in mst.values()):
        fam["notes"].append(
            "--b1-unpol mstw (the configuration the A = 2 gate passes on): A_zz = "
            "%.7e (unset) / %.7e (r1998), ratio %.6f; sigma %.10g / %.10g pb." % (
                mst["unset"]["azz"], mst["r1998"]["azz"],
                mst["r1998"]["azz"] / mst["unset"]["azz"],
                mst["unset"]["sigma"], mst["r1998"]["sigma"]))
        fam["b1_unpol_mstw"] = dict(azz_unset=mst["unset"]["azz"],
                                    azz_r1998=mst["r1998"]["azz"],
                                    sigma_unset=mst["unset"]["sigma"],
                                    sigma_r1998=mst["r1998"]["sigma"])
    else:
        fam["notes"].append("--b1-unpol mstw: not run -- %s" % "; ".join(
            str(v) for v in mst.values() if not isinstance(v, dict)))
    # the shipped miller b1 at the same window: does b1_model move A2?
    out, err = _attempt(dict(isotope="6Li", channel="inclusive", scenario=sc,
                             plan="transverse-tensor"), events, seed)
    if err is None:
        v, _e, _c = _superbin_amplitudes(out)
        fam["miller_A2"] = [v["A2_S%d" % i] for i in range(1, 5)]
        fam["notes"].append(
            "b1_model at the same window (x_max %.2f): A2(S1..S4) li6-convolution / "
            "miller = %s -- the b1 term reaches A2 only through the cell weights "
            "(1 + w_avg)." % (X_MAX_CONVOLUTION, " / ".join(
                "%.6f" % (default["values"]["A2_S%d" % i] / fam["miller_A2"][i - 1])
                for i in range(1, 5))))
        del out
    # the documented function-level points (README.md / docs/USAGE.md):
    # x = 0.05 / 0.10 / 0.30 at y = 0.5, Q2 = 2.5 -- not a point of config 1
    if "unset" in kernels and "r1998" in kernels:
        K = lg.InclusiveKernel
        pts = []
        for x in (0.05, 0.10, 0.30):
            q2, y = 2.5, 0.5
            res = {}
            for key, k in kernels.items():
                t = k.tables(x, q2)
                d = K.dphi(t, x, y)
                res[key] = (K.tensor_kernel(t, x, y) / d, t.delta / d)
            pts.append(dict(x=x, kd=100.0 * (res["r1998"][0] / res["unset"][0] - 1.0),
                            amp=100.0 * (res["r1998"][1] / res["unset"][1] - 1.0)))
        fam["function_points"] = pts
        # why the sweet spots move by ~ -15 %: at their small y, D_phi is
        # dominated by (1-y)/(x y^2) F2, so A2 ~ x Delta/F2 ~ F1/F2 ~ 1/(1 + R)
        rr = [(sp[0], float(lg.r_sigma_lt(sp[1], sp[2])), float(lg.r1998(sp[1], sp[2])))
              for sp in SWEET_SPOTS]
        fam["r_at_centres"] = rr
        fam["notes"].append(
            "R at the super-bin centres: r_sigma_lt %s, r1998 %s; (1 + R_sigma_lt) / "
            "(1 + R_1998) = %s -- the small-y limit of the A2 ratio (D_phi = F1 + "
            "(1-y)/(x y^2) F2 is F2-dominated there, so A2 ~ F1/F2 ~ 1/(1 + R))." % (
                " / ".join("%.4f" % r[1] for r in rr), " / ".join("%.4f" % r[2] for r in rr),
                " / ".join("%.4f" % ((1.0 + r[1]) / (1.0 + r[2])) for r in rr)))
        fam["notes"].append(
            "the documented FUNCTION-LEVEL points (README.md --r-source row, "
            "docs/USAGE.md): at x = 0.05 / 0.10 / 0.30, y = 0.5, Q2 = 2.5 the shared "
            "r1998 moves K/D_phi by %s %% and Delta/D_phi (the cos 2phi amplitude) "
            "by %s %% -- reproduced; those points are not the sweet spots." % (
                " / ".join("%+.4f" % p["kd"] for p in pts),
                " / ".join("%+.4f" % p["amp"] for p in pts)))
    return fam


def family_li7_inclusive(events, seed, log):
    fam = _family(
        "li7-incl",
        "Inclusive 7Li -- A_par of paper 3 (rank 5)",
        "7Li, channel inclusive, plan helicity-flip (P_e = 0.7, P_z = 0.7, "
        "ladder fill); x_max 1.0; N = %d, seed %d" % (events, seed),
        [("sigma", "sigma [pb]", "pb"), ("apar", "A_par", "asym")])
    base = dict(isotope="7Li", channel="inclusive", plan="helicity-flip", pe=PE)

    def ext(out):
        a = _apar(out)
        return dict(sigma=a["sigma"], apar=a["apar"]), {}

    for pol in ("toy", "nnpdfpol"):
        for unpol in ("toy", "mstw", "ct18nlo"):
            kw = {}
            if pol != "toy":
                kw["pol_sf"] = pol
            if unpol != "toy":
                kw["unpol_sf"] = unpol
            setting = "pol_sf = %s, unpol_sf = %s%s" % (
                pol, unpol, " (shipped default)" if not kw else "")
            knob = None if not kw else ("pol_sf" if pol != "toy" else "unpol_sf")
            out = _price(fam, setting, knob, "%s/%s" % (pol, unpol), dict(base, **kw),
                         ext, events, seed, log)
            if out is not None and pol != "toy" and unpol != "toy":
                if _kp(out, "unpol_sf")["status"] != "read":
                    raise RuntimeError("unpol_sf not read on %s" % setting)
            del out
    default = fam["rows"][0]
    for knob, val, expect in (("r_source", "r1998", "refused"),
                              ("cluster_wave", "vmc", "refused"),
                              ("triton_sf", "ciofi-simula", "not-read"),
                              ("eps_b0", EPS_B0_MEASURED, "not-read"),
                              ("f0", F0_BAND[0], "not-read")):
        fam["rows"].append(_probe(base, knob, val, events, seed, default, ext, expect))
    _ratios(fam)
    return fam


def family_li7_tagged(events, seed, log):
    fam = _family(
        "li7-tag",
        "Tagged 7Li alpha -- A_par on the quasi-free triton and the spectator "
        "<k> (paper 2, rank 3)",
        "7Li, channel tagged-7Li-alpha, plan helicity-flip (P_e = 0.7, P_z = 0.7, "
        "ladder fill); optics yr-high-acceptance; N = %d, seed %d" % (events, seed),
        [("sigma", "sigma [pb]", "pb"), ("apar", "A_par", "asym"),
         ("k_mean", "<k> generated [GeV]", "k"),
         ("rp_frac", "RP tag fraction (generated)", "frac")])
    base = dict(isotope="7Li", channel="tagged-alpha", plan="helicity-flip", pe=PE)
    # the T1 final-state columns a READ --triton-sf moves, kept to say WHAT it moves
    final_state = ("n_partner", "struck_virtuality", "struck_pol")
    kept = {}
    out = _price(fam, "shipped default (%s, hulthen, toy, toy)" % _waves()[0], None,
                 None, base, _tagged_li7, events, seed, log)
    if out is not None:
        kept["default"] = {c: np.asarray(out[c]).copy() for c in final_state}
    del out
    for knob, val in (("cluster_wave", _waves()[1]), ("triton_sf", "ciofi-simula"),
                      ("pol_sf", "nnpdfpol"), ("unpol_sf", "mstw"),
                      ("unpol_sf", "ct18nlo")):
        out = _price(fam, "%s = %s" % (knob, val), knob, val,
                     dict(base, **{knob: val}), _tagged_li7, events, seed, log)
        if out is not None and knob == "triton_sf":
            kept["triton"] = {c: np.asarray(out[c]).copy() for c in final_state}
        del out
    if "default" in kept and "triton" in kept:
        d, t = kept["default"], kept["triton"]
        fam["triton_final_state"] = {
            c: dict(n_differ=int(np.sum(d[c] != t[c])), mean_default=float(d[c].mean()),
                    mean_ciofi=float(t[c].mean())) for c in final_state}
        fam["notes"].append(
            "--triton-sf is READ here, and on the same seed it moves the T1 final "
            "state, not the ranked observables: %s; the per-category cross sections "
            "(so A_par) are bit-identical and <k> moves by its sampling only." % "; ".join(
                "%s differs on %d of %d events (mean %.6g -> %.6g)" % (
                    c, v["n_differ"], events, v["mean_default"], v["mean_ciofi"])
                for c, v in fam["triton_final_state"].items()))
    default = fam["rows"][0]
    for knob, val, expect in (("r_source", "r1998", "refused"),
                              ("eps_b0", EPS_B0_MEASURED, "not-read"),
                              ("f0", F0_BAND[0], "not-read")):
        fam["rows"].append(_probe(base, knob, val, events, seed, default,
                                  _tagged_li7, expect))
    _ratios(fam)
    return fam


def family_li6_tagged(events, seed, log):
    fam = _family(
        "li6-tag",
        "Tagged 6Li alpha -- A_zz^tag(k) on the embedded deuteron (paper 2, rank 6)",
        "6Li, channel tagged-6Li-alpha, plan tensor-thirds (P_z = 0.7, P_zz = 0.6); "
        "A_zz^tag = azz_tensor_curve_weighted at acceptance_weights(yr-high-acceptance, "
        "n_phi = %d), at the TaggedModel cell nearest k; N = %d, seed %d"
        % (N_PHI_ACCEPTANCE, events, seed),
        [("sigma", "sigma [pb]", "pb"),
         ("azz_tag_020", "A_zz^tag(k ~ 0.20)", "azz"),
         ("azz_tag_030", "A_zz^tag(k ~ 0.30)", "azz"),
         ("rp_frac", "RP tag fraction (generated)", "frac")])
    base = dict(isotope="6Li", channel="tagged-alpha", plan="tensor-thirds")
    _price(fam, "shipped default (cluster_wave = %s)" % _waves()[0], None, None, base,
           _tagged_li6, events, seed, log)
    for knob, val in (("cluster_wave", _waves()[1]), ("unpol_sf", "mstw"),
                      ("unpol_sf", "ct18nlo")):
        _price(fam, "%s = %s" % (knob, val), knob, val, dict(base, **{knob: val}),
               _tagged_li6, events, seed, log)
    default, (dw, aw) = fam["rows"][0], _waves()
    alt = [r for r in fam["rows"] if r["setting"] == "cluster_wave = %s" % aw
           and r["provenance"] == "read"]
    if alt:
        a, b = default["values"], alt[0]["values"]
        fam["notes"].append(
            "the k- and angle-integrated struck-deuteron dilutions (TaggedModel."
            "vector_dilution / tensor_dilution, N-independent): %s %.6f / %.6f, "
            "%s %.6f / %.6f (x%.5f / x%.5f).  A_zz on sigma_tot is zero by "
            "construction on this channel (%.1e / %.1e), and with NO acceptance "
            "weight the same curve reads %.1e / %.1e (%s) and %.1e / %.1e (%s) "
            "at the two k cells: the tensor signal is the cos theta_k shape of the "
            "momentum distribution, so A_zz^tag(k) exists only through the "
            "acceptance." % (
                dw, a["vec_dil"], a["tens_dil"], aw, b["vec_dil"], b["tens_dil"],
                b["vec_dil"] / a["vec_dil"], b["tens_dil"] / a["tens_dil"],
                a["azz_sigma_tot"], b["azz_sigma_tot"], a["azz_4pi_020"],
                a["azz_4pi_030"], dw, b["azz_4pi_020"], b["azz_4pi_030"], aw))
    for knob, val, expect in (("pol_sf", "nnpdfpol", "not-read"),
                              ("triton_sf", "ciofi-simula", "not-read"),
                              ("r_source", "r1998", "refused"),
                              ("eps_b0", EPS_B0_MEASURED, "not-read"),
                              ("f0", F0_BAND[0], "not-read")):
        fam["rows"].append(_probe(base, knob, val, events, seed, default,
                                  _tagged_li6, expect))
    _ratios(fam)
    return fam


def family_li6_coherent(events, seed, log):
    fam = _family(
        "li6-coh",
        "Coherent 6Li -- c2 of 1 + c2 cos 2(phi_t - phi_S) (claim 2, rank 2)",
        "6Li, channel coherent, plan transverse-tensor (P_zz = 0.6); scenario "
        "defaults f0 = 0.04, slope_b = 50, amp = 0.01, eps_b0 = -0.08 unless the "
        "row says; |t| <= 0.2 GeV^2; N = %d, seed %d" % (events, seed),
        [("sigma", "sigma_coh [pb]", "pb"),
         ("c2_exp", "c2 at the spectrum's mean abs(t) (exact)", "c2"),
         ("c2_tag", "c2 at the tag mean abs(t) = 0.06 GeV^2", "c2"),
         ("c2_gen", "c2 generated (sampled; no ratio)", "c2")])
    base = dict(isotope="6Li", channel="coherent", plan="transverse-tensor")
    _price(fam, "shipped default (eps_b0 = -0.08, f0 = 0.04)", None, None, base,
           _coherent, events, seed, log)
    for knob, val in (("eps_b0", EPS_B0_MEASURED), ("f0", F0_BAND[0]),
                      ("f0", F0_BAND[1])):
        _price(fam, "%s = %g (API: coherent={%r: %g})" % (knob, val, knob, val)
               if knob == "eps_b0" else "%s = %g (--coherent-f0)" % (knob, val),
               knob, val, dict(base, coherent={knob: val}), _coherent, events, seed, log)
    for val in ("mstw", "ct18nlo"):
        _price(fam, "unpol_sf = %s" % val, "unpol_sf", val, dict(base, unpol_sf=val),
               _coherent, events, seed, log)
    default = fam["rows"][0]
    dv = default["values"]
    fam["notes"].append(
        "c2 = P_zz (-eps_b0 B |t| / 2 + amp) is linear in |t| (coherent.hpp), so its "
        "expectation over the generated spectrum is c2 at the spectrum's exact <|t|> "
        "= %.6f GeV^2 (B exp(-B|t|) on [0, %.2f], B = %.0f); the tag <|t|> = %.4f "
        "GeV^2 is CoherentScenario::mean_t_tagged(%.2f GeV), the letter's envelope; "
        "%d of N coherent recoils reach the pots at the YR envelope." % (
            dv["t_exp"], dv["t_max"], dv["slope_b"], dv["t_tag"], TAG_ENVELOPE_PT,
            dv["rp_n"]))
    for knob, val, expect in (("pol_sf", "nnpdfpol", "not-read"),
                              ("cluster_wave", "vmc", "refused"),
                              ("triton_sf", "ciofi-simula", "not-read"),
                              ("r_source", "r1998", "refused")):
        fam["rows"].append(_probe(base, knob, val, events, seed, default,
                                  _coherent, expect))
    _ratios(fam)
    return fam


BUILDERS = {"li6-incl": family_li6_inclusive, "li6-incl-r": family_li6_inclusive_r,
            "li7-incl": family_li7_inclusive, "li7-tag": family_li7_tagged,
            "li6-tag": family_li6_tagged, "li6-coh": family_li6_coherent}


def price_table(events=N_DEFAULT, seed=SEED, families=FAMILIES, log=None):
    """Build the table: {'events', 'seed', 'families': {name: family}}."""
    log = log or (lambda s: None)
    t0 = time.time()
    fams = {}
    for name in families:
        if name not in BUILDERS:
            raise ValueError("unknown family %r; know %s" % (name, ", ".join(FAMILIES)))
        t = time.time()
        fams[name] = BUILDERS[name](events, seed, log)
        fams[name]["seconds"] = round(time.time() - t, 2)
    return dict(name=NAME, events=events, seed=seed, families=fams,
                lipolgen=lg.__version__, seconds=round(time.time() - t0, 2),
                env=dict(lhapdf=bool(lg.HAVE_LHAPDF), pythia=bool(lg.HAVE_PYTHIA8)))


# ------------------------------------------------------------- printing

def _fmt(v, kind):
    if v is None or (isinstance(v, float) and not math.isfinite(v)):
        return "--"
    if kind == "pb":
        return "%.10g" % v
    if kind in ("asym",):
        return "%.6e" % v
    if kind == "azz":
        return "%.4f" % v
    if kind == "c2":
        return "%.6f" % v
    if kind == "k":
        return "%.5f" % v
    if kind == "frac":
        return "%.5f" % v
    return "%.6g" % v


def _cell(row, key, kind, default):
    v = row["values"].get(key)
    if row["provenance"] not in ("read", "default"):
        return "--"
    s = _fmt(v, kind)
    e = row["errors"].get(key)
    if e is not None and math.isfinite(e):
        s += " +- " + _fmt(e, kind)
    r = row["ratios"].get(key)
    if r is not None and not default and key not in NO_RATIO:
        s += " (x%.5f)" % r
    return s


def markdown(table):
    lines = ["# %s -- N = %d, seed %d, lipolgen %s" % (
        NAME, table["events"], table["seed"], table["lipolgen"]), ""]
    foot = []
    for fam in table["families"].values():
        lines += ["## " + fam["title"], "", "*Window:* " + fam["window"] + ".", ""]
        lines.append("| setting | provenance | " + " | ".join(c[1] for c in fam["columns"]) + " |")
        lines.append("|---|---|" + "---|" * len(fam["columns"]))
        for i, row in enumerate(fam["rows"]):
            prov = row["provenance"]
            if prov == "not-read":
                prov = "not-read: %s; bit-identical to the default: %s" % (
                    row["label"], "yes" if row["identical"] else "NO")
            elif prov in ("refused", "environment"):
                text = " ".join(row["label"].split())
                if text not in foot:
                    foot.append(text)
                prov = "%s [%d]" % (prov, foot.index(text) + 1)
            cells = [_cell(row, c[0], c[2], row["provenance"] == "default")
                     for c in fam["columns"]]
            lines.append("| %s | %s | %s |" % (row["setting"], prov, " | ".join(cells)))
        lines.append("")
        for n in fam["notes"]:
            lines.append("- " + n)
        if fam["notes"]:
            lines.append("")
    if foot:
        lines.append("Refusal and environment texts, verbatim:")
        lines.append("")
        for i, f in enumerate(foot, 1):
            lines.append("[%d] %s" % (i, f))
        lines.append("")
    return "\n".join(lines)


def _closure_lines(table):
    """The generated columns against the exact numbers (the default rows)."""
    out = []
    f = table["families"]
    if "li6-incl" in f:
        d = f["li6-incl"]["rows"][0]
        for i in range(1, 5):
            v = d["values"]
            if "gen_A2_S%d" % i in v:
                e = d["errors"]["gen_A2_S%d" % i]
                out.append("li6-incl S%d: exact A2 %.6e; generated moment on %d events "
                           "%.4e +- %.4e (pull %+.2f)" % (
                               i, v["A2_S%d" % i], v["n_S%d" % i], v["gen_A2_S%d" % i], e,
                               (v["gen_A2_S%d" % i] - v["A2_S%d" % i]) / e))
    if "li6-tag" in f:
        for r in f["li6-tag"]["rows"][:2]:
            if r["provenance"] not in ("read", "default"):
                continue
            for tag in ("020", "030"):
                v, e = r["values"], r["errors"]
                if "gen_azz_" + tag in v:
                    out.append("li6-tag %s k bin %s +- %.3f: model bin average %.4f; "
                               "generated RP-accepted n = %s, A_zz %.4f +- %.4f (pull %+.2f)" % (
                                   r["setting"], "0.%s" % tag[1:], K_HALF_WIDTH,
                                   v["azz_tag_bin_" + tag], v["gen_n_" + tag],
                                   v["gen_azz_" + tag], e["gen_azz_" + tag],
                                   (v["gen_azz_" + tag] - v["azz_tag_bin_" + tag]) / e["gen_azz_" + tag]))
    if "li6-coh" in f:
        for r in f["li6-coh"]["rows"][:2]:
            if r["provenance"] not in ("read", "default"):
                continue
            v, e = r["values"], r["errors"]
            out.append("li6-coh %s: generated c2 %.4f +- %.4f against the model %.5f over "
                       "the same |t| sample (pull %+.2f) and the exact %.5f; generated <|t|> "
                       "%.6f +- %.6f against the exact %.6f GeV^2; RP-accepted %d of N" % (
                           r["setting"], v["c2_gen"], e["c2_gen"], v["c2_sample"],
                           (v["c2_gen"] - v["c2_sample"]) / e["c2_gen"], v["c2_exp"],
                           v["t_mean"], e["t_mean"], v["t_exp"], v["rp_n"]))
    return out


# ------------------------------------------------------------- the figure

def _r(table, fam, setting_prefix, key):
    f = table["families"].get(fam)
    if not f:
        return None
    for row in f["rows"]:
        if row["setting"].startswith(setting_prefix) and row["provenance"] == "read":
            return row["ratios"].get(key)
    return None


def plot(table, out_dir=None):
    """docs/benchmarking/figures/price_defaults.png (or out_dir/): the ratio
    to the shipped default, per knob, in the fixed slot colours of
    validation/benchmarks/_plots.py; top the nucleon structure-function
    inputs (ingredients A and H), bottom the nuclear and coherent-scenario
    inputs (G, E).  Only priced (READ) cells are drawn; not-read, refused and
    environment cells are the table's, never a point."""
    if BENCH not in sys.path:
        sys.path.insert(0, BENCH)
    import _plots
    plt = _plots.mpl()
    fig, (ax_a, ax_b) = plt.subplots(2, 1, figsize=(6.4, 7.2), layout="constrained")
    fig.get_layout_engine().set(rect=(0.0, 0.03, 1.0, 0.905), hspace=0.06)

    cats_a = [("li6-incl", "sigma", r"$^6$Li $\sigma$"),
              ("li6-incl", "A2_S1", r"$A_2$ S1"), ("li6-incl", "A2_S2", r"$A_2$ S2"),
              ("li6-incl", "A2_S3", r"$A_2$ S3"), ("li6-incl", "A2_S4", r"$A_2$ S4"),
              ("li6-incl", "azz", r"$^6$Li $A_{zz}$"),
              ("li7-incl", "sigma", r"$^7$Li $\sigma$"),
              ("li7-incl", "apar", r"$^7$Li $A_\parallel$"),
              ("li7-tag", "apar", r"$^7$Li-$\alpha$ $A_\parallel$"),
              ("li6-coh", "sigma", r"coh $\sigma$"),
              ("li6-coh", "c2_exp", r"coh $c_2$")]
    series_a = [("mstw", "unpol_sf = mstw", "o", "unpol_sf = mstw"),
                ("ct18nlo", "unpol_sf = ct18nlo", "s", "unpol_sf = ct18nlo"),
                ("nnpdfpol", "pol_sf = nnpdfpol", "D", None),
                ("r1998", "r_source = r1998 (li6-convolution path)", "^", None)]

    def ratio_a(fam, key, slot):
        if slot == "mstw":
            return _r(table, fam, "unpol_sf = mstw" if fam != "li7-incl"
                      else "pol_sf = toy, unpol_sf = mstw", key)
        if slot == "ct18nlo":
            return _r(table, fam, "unpol_sf = ct18nlo" if fam != "li7-incl"
                      else "pol_sf = toy, unpol_sf = ct18nlo", key)
        if slot == "nnpdfpol":
            return _r(table, fam, "pol_sf = nnpdfpol" if fam != "li7-incl"
                      else "pol_sf = nnpdfpol, unpol_sf = toy", key)
        if slot == "r1998":
            if fam != "li6-incl":
                return None
            return _r(table, "li6-incl-r", "r_source = r1998", key)
        return None

    cats_b = [("li7-tag", "apar", r"$^7$Li-$\alpha$ $A_\parallel$"),
              ("li7-tag", "k_mean", r"$^7$Li-$\alpha$ $\langle k\rangle$"),
              ("li7-tag", "rp_frac", r"$^7$Li-$\alpha$ tag frac."),
              ("li6-tag", "azz_tag_020", r"$A_{zz}^{tag}$(0.20)"),
              ("li6-tag", "azz_tag_030", r"$A_{zz}^{tag}$(0.30)"),
              ("li6-tag", "rp_frac", r"$^6$Li-$\alpha$ tag frac."),
              ("li6-coh", "sigma", r"coh $\sigma$"),
              ("li6-coh", "c2_exp", r"coh $c_2$ ($\langle|t|\rangle$ = 0.020)"),
              ("li6-coh", "c2_tag", r"coh $c_2$ (tag $\langle|t|\rangle$)")]
    series_b = [("alt", "cluster_wave = " + _waves()[1], "o", "cluster_wave = " + _waves()[1]),
                ("alt2", "triton_sf = ciofi-simula", "s", "triton_sf = ciofi-simula"),
                ("alt3", "eps_b0 = -0.0070 (measured Q)", "D", "eps_b0 = ")]

    def draw(ax, cats, series, getter, f0_band=False):
        xs = np.arange(len(cats))
        offs = np.linspace(-0.24, 0.24, len(series))
        drawn = set()
        drawn_y = []
        for (slot, label, marker, _p), off in zip(series, offs):
            px, py = [], []
            for i, (fam, key, _lab) in enumerate(cats):
                r = getter(fam, key, slot)
                if r is not None and r > 0:
                    px.append(xs[i] + off)
                    py.append(r)
            if px:
                ax.plot(px, py, **_plots.style_series(slot, label=label, marker=marker))
                drawn.add(slot)
                drawn_y.extend(py)
        if f0_band:
            lo = _r(table, "li6-coh", "f0 = %g" % F0_BAND[0], "sigma")
            hi = _r(table, "li6-coh", "f0 = %g" % F0_BAND[1], "sigma")
            c2lo = _r(table, "li6-coh", "f0 = %g" % F0_BAND[0], "c2_exp")
            c2hi = _r(table, "li6-coh", "f0 = %g" % F0_BAND[1], "c2_exp")
            first = True
            for i, (fam, key, _l) in enumerate(cats):
                if fam != "li6-coh":
                    continue
                pair = (lo, hi) if key == "sigma" else (c2lo, c2hi) if key == "c2_exp" else (None, None)
                if None in pair:
                    continue
                ax.fill_between([xs[i] - 0.34, xs[i] - 0.14], min(pair), max(pair),
                                color=_plots.BAND, alpha=0.9, linewidth=0, zorder=1,
                                label="f0 in {0.02, 0.08} (band)" if first else None)
                drawn_y.extend(pair)
                first = False
        ax.axhline(1.0, color=_plots.INK, linewidth=0.8, zorder=2)
        ax.set_yscale("log")
        # F6 (2026-10-29): the range is the drawn points' own, never narrower
        # than 0.1 .. 2.6.  It was FIXED at 0.1 .. 2.6, so once D6 made hulthen
        # the priced cluster wave its A_zz^tag(0.30) ratio, x4.78, had no marker.
        lo = min([0.1] + [v / 1.25 for v in drawn_y])
        hi = max([2.6] + [v * 1.25 for v in drawn_y])
        ticks = [t for t in (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0, 10.0,
                             20.0, 50.0) if lo <= t <= hi]
        ax.set_yticks(ticks)
        ax.set_yticklabels(["%g" % t for t in ticks])
        ax.set_ylim(lo, hi)
        ax.set_xticks(xs)
        ax.set_xticklabels([c[2] for c in cats], rotation=30, ha="right", fontsize=7)
        ax.set_ylabel("ratio to the shipped default")
        ax.legend(loc="lower left", ncol=2, fontsize=7)

    draw(ax_a, cats_a, series_a, ratio_a)
    ax_a.set_title("nucleon structure-function inputs (ingredients A, H)",
                   fontsize=8, loc="left", color=_plots.INK_SECONDARY)

    def ratio_b(fam, key, slot):
        if slot == "alt":
            return _r(table, fam, "cluster_wave = " + _waves()[1], key)
        if slot == "alt2":
            return _r(table, fam, "triton_sf = ciofi-simula", key)
        if slot == "alt3":
            return _r(table, fam, "eps_b0 = ", key)
        return None

    draw(ax_b, cats_b, series_b, ratio_b, f0_band=True)
    ax_b.set_title("nuclear-model and coherent-scenario inputs (ingredients G, E)",
                   fontsize=8, loc="left", color=_plots.INK_SECONDARY)
    _plots.status_title(
        fig, NAME, "priced",
        "config 1 (e 10 x ion 99.5 GeV/u), N = %d, seed %d; one knob at a time, "
        "same channel and plan" % (table["events"], table["seed"]))
    return _plots.save(fig, NAME, out_dir)


# ----------------------------------------------------------------- main

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--events", type=int, default=N_DEFAULT)
    ap.add_argument("--seed", type=int, default=SEED)
    ap.add_argument("--families", default=",".join(FAMILIES),
                    help="comma-separated subset of %s" % ", ".join(FAMILIES))
    ap.add_argument("--plot", action="store_true")
    ap.add_argument("--out-dir", default=None,
                    help="where --plot writes price_defaults.png "
                         "(default docs/benchmarking/figures/)")
    ap.add_argument("--json", default=None, help="also write the table as JSON here")
    args = ap.parse_args(argv)
    fams = tuple(f for f in args.families.split(",") if f)
    try:
        table = price_table(args.events, args.seed, fams,
                            log=lambda s: print(s, file=sys.stderr))
        print(markdown(table))
        cl = _closure_lines(table)
        if cl:
            print("Closure of the exact numbers against the generated columns:")
            print("")
            for c in cl:
                print("- " + c)
            print("")
        print("wall time %.1f s (%s)" % (table["seconds"], ", ".join(
            "%s %.1f s" % (k, v["seconds"]) for k, v in table["families"].items())))
        if args.json:
            with open(args.json, "w") as fh:
                json.dump(table, fh, indent=1, default=float)
            print("wrote %s" % args.json)
        if args.plot:
            print("wrote %s" % plot(table, args.out_dir))
    except Exception:  # the harness convention: 2 = broken
        import traceback
        traceback.print_exc()
        return 2
    return 0


def _waves():
    """D6 (2026-10-27): (the shipped --cluster-wave, the other one), READ from
    the library's `CLUSTER_WAVE_DEFAULT` (cluster.hpp) -- so the tagged
    families' default row is named by what ships, and the priced alternative
    is the wave that does not (hulthen since D6; vmc before)."""
    d = lg._lipolgen.cluster_wave_name(lg._lipolgen.CLUSTER_WAVE_DEFAULT)
    return d, ("hulthen" if d == "vmc" else "vmc")


if __name__ == "__main__":
    sys.exit(main())
