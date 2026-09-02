#!/usr/bin/env python3
"""What switching the cluster radial forms to the ANL VMC tables does to the
TAGGED observables: Roman-Pot tag fractions and the tagged tensor asymmetry.

Run after `validation/vmc_reconcile.py --markdown docs/open_items/vmc_reconciliation.md`;
with `--markdown` it appends its tables to that file after the marker the
reconciliation script leaves.

Two things are measured, both for 6Li α+d and 7Li α+t, VMC vs Hulthen
β = 0.30 (the model default):

1. **Tag fraction** — the fraction of generated tagged events whose spectator
   lands in the Roman Pots (`export.rp_accepted`: main window + near-beam
   tail), at BOTH the Yellow Report high-acceptance envelope and the tagging
   optics, e × ion = 10 × 99.5 GeV/u (`--config 1`).  This is exactly the
   number `python -m lipolgen.cli` prints as "tag fraction at <optics>".

2. **A_zz^tag(k)** — `azz_tensor_curve_weighted`, the acceptance-weighted
   wave-function tensor asymmetry (n_+1 + n_-1 - 2 n_0)/(n_+1 + n_-1 + n_0)
   versus the spectator momentum k.  6Li only: it needs a spin-1 channel.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
MARKER = "<!-- APPEND: validation/vmc_tag_fractions.py -->"

OPTICS = ("yr-high-acceptance", "tagging")
WAVES = ("hulthen", "vmc")
K_ROWS = (0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.40, 0.50)


def load():
    build = ROOT / "build" / "python"
    if build.is_dir():
        sys.path.insert(0, str(build))
    import lipolgen
    return lipolgen


def tag_fractions(lp, isotope, n_events, seed, config=1):
    """{(optics, wave): (all, main, near-beam)} tag fractions."""
    from lipolgen import export
    out = {}
    for optics in OPTICS:
        for wave in WAVES:
            cfg = lp.make_config(isotope=isotope, config=config,
                                 channel="tagged-alpha", events=n_events,
                                 seed=seed, optics=optics, cluster_wave=wave)
            # tensor-thirds is spin-1 only; 7Li (J = 3/2) needs helicity-flip
            j = lp.ion_spin(cfg.isotope)
            plan = lp.make_plan("tensor-thirds" if j == 1.0 else "helicity-flip",
                                j=j)
            p = lp._lipolgen.Pipeline(cfg, plan)
            cols = p.generate(0, n_events, 1)
            route = cols["route"]
            out[(optics, wave)] = (
                float(np.mean(export.rp_accepted(cols))),
                float(np.mean(route == lp._lipolgen.Route.RomanPots)),
                float(np.mean(route == lp._lipolgen.Route.RPNearBeam)),
                p.optics.name,
            )
    return out


def beam_energies(lp, config):
    """{isotope: (electron GeV, ion GeV/u)} at one `default_configs` index."""
    L = lp._lipolgen
    return {iso: (c.electron_energy, c.ion_momentum_per_nucleon)
            for iso in ("6Li", "7Li")
            for c in (L.default_configs(iso)[config],)}


def energy_label(energies):
    """'e x ion = 10 x 99.5 GeV/u', or a per-isotope breakdown when the
    rigidity cap makes the two isotopes' ion energy differ (config 2)."""
    vals = set(energies.values())
    if len(vals) == 1:
        e, u = next(iter(vals))
        return f"e x ion = {e:g} x {u:g} GeV/u"
    return "e x ion: " + ", ".join(
        f"{iso} {e:g} x {u:g} GeV/u" for iso, (e, u) in energies.items())


def li6_variants(lp):
    """The four 6Li channels the A_zz comparison needs, in report order.

    `hulthen`         the model default, beta = 0.30, P_D = 0.0867
    `hulthen-P_D_vmc` the same analytic shapes at the VMC D-state probability
                      -- isolates the P_D change from the SHAPE change
    `vmc`             the production VMC channel
    `vmc-flipD`       the VMC channel with the D-wave table negated -- the
                      control on the one input that the momentum files alone
                      cannot supply, the S-D RELATIVE SIGN
    """
    L = lp._lipolgen
    out = {}
    out["hulthen"] = L.li6_alpha_channel()
    out["hulthen-P_D_vmc"] = L.li6_alpha_channel(p_d=L.VMC_P_D_LI6)
    vmc = L.li6_alpha_channel(source=L.ClusterWaveSource.VmcAV18)
    out["vmc"] = vmc
    flipped = L.li6_alpha_channel(source=L.ClusterWaveSource.VmcAV18)
    waves = list(flipped.waves)
    d = waves[1]
    d.vmc = L.VmcRadial(list(np.array(d.vmc.k)),
                        list(-np.array(d.vmc.psi)), d.l,
                        d.vmc.provenance + " [D-WAVE SIGN FLIPPED: control]")
    waves[1] = d
    flipped.waves = waves
    out["vmc-flipD"] = flipped
    return out


def azz_curves(lp):
    """{variant: (k, A_zz^tag(k))} for 6Li at the YR high-acceptance optics."""
    L = lp._lipolgen
    beams = L.default_configs("6Li")[1]
    p_u = beams.ion_momentum_per_nucleon
    optics = L.yr_optics("6Li", p_u, True)
    pot = L.yr_config_key("6Li", p_u)
    out = {}
    for name, ch in li6_variants(lp).items():
        m = L.TaggedModel(ch)
        eps = L.acceptance_weights(m, p_u, optics, pot, 64)
        out[name] = (np.array(m.k), np.array(L.azz_tensor_curve_weighted(m, eps)))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--events", type=int, default=60000)
    ap.add_argument("--seed", type=int, default=20260829)
    ap.add_argument("--markdown", type=Path, default=None)
    ap.add_argument("--configs", type=str, default="1",
                    help="comma-separated default_configs() indices for the "
                         "Roman-Pot tag-fraction section, e.g. 0,1,2 "
                         "(default: 1, the 10 x 99.5 GeV/u point)")
    args = ap.parse_args()
    configs = [int(c) for c in args.configs.split(",")]

    lp = load()
    L = lp._lipolgen

    # ---- what actually changed in the model -----------------------------
    print("=" * 92)
    print("0.  THE CHANNELS")
    print("=" * 92)
    for iso, maker in (("6Li", L.li6_alpha_channel), ("7Li", L.li7_alpha_channel)):
        for wave in WAVES:
            src = (L.ClusterWaveSource.Hulthen if wave == "hulthen"
                   else L.ClusterWaveSource.VmcAV18)
            ch = maker(source=src) if iso == "7Li" else maker(source=src)
            print(f"{iso} {wave:<8} {ch.label}")
            for w in ch.waves:
                extra = "" if w.vmc is None else ("  <- " + w.vmc.provenance)
                print(f"      L={w.l}  P_L={w.prob:.6f}  beta={w.beta:.2f}{extra}")

    # ---- moments on the TaggedModel grid --------------------------------
    print()
    print("=" * 92)
    print("1.  MOMENTS ON THE TAGGED-MODEL GRID (k in [1e-4, 1.2], 280 cells)")
    print("=" * 92)
    print(f"{'channel / waves':<30}{'<k> GeV':>10}{'P(k>0.2)':>11}"
          f"{'P(k>0.3)':>11}{'P(k>0.45)':>11}{'P_D':>10}")
    moments = {}
    for iso, maker in (("6Li", L.li6_alpha_channel), ("7Li", L.li7_alpha_channel)):
        for wave in WAVES:
            src = (L.ClusterWaveSource.Hulthen if wave == "hulthen"
                   else L.ClusterWaveSource.VmcAV18)
            ch = maker(source=src)
            m = L.TaggedModel(ch)
            k = np.array(m.k)
            dens = np.zeros_like(k)
            d2 = np.zeros_like(k)
            for w in ch.waves:
                r = np.array(m.radial_table(w.l))
                dens = dens + k ** 2 * r ** 2
                if w.l == 2:
                    d2 = d2 + k ** 2 * r ** 2
            tot = np.trapz(dens, k)
            dens = dens / tot
            mk = np.trapz(k * dens, k)

            def tail(k0):
                msk = k >= k0
                return float(np.trapz(dens[msk], k[msk]))

            pd = float(np.trapz(d2, k) / tot) if d2.any() else float("nan")
            moments[(iso, wave)] = (mk, tail(0.2), tail(0.3), tail(0.45), pd)
            pd_s = "     --" if np.isnan(pd) else f"{pd:10.5f}"
            print(f"{iso + ' ' + wave:<30}{mk:10.4f}{tail(0.2):11.4f}"
                  f"{tail(0.3):11.4f}{tail(0.45):11.4f}{pd_s}")

    # ---- tag fractions ---------------------------------------------------
    tf_by_cfg = {}
    for ci, cfg_idx in enumerate(configs):
        label = energy_label(beam_energies(lp, cfg_idx))
        header = "2" if len(configs) == 1 else f"2{chr(ord('a') + ci)}"
        print()
        print("=" * 92)
        print(f"{header}.  ROMAN-POT TAG FRACTIONS ({args.events} events, "
              f"{label}, seed {args.seed})")
        print("=" * 92)
        tf = {}
        for iso in ("6Li", "7Li"):
            tf[iso] = tag_fractions(lp, iso, args.events, args.seed,
                                    config=cfg_idx)
        tf_by_cfg[cfg_idx] = tf
        print(f"{'isotope':<8}{'optics':<22}{'envelope':<26}"
              f"{'Hulthen b=0.30':>16}{'VMC AV18':>12}{'ratio':>9}")
        for iso in ("6Li", "7Li"):
            for optics in OPTICS:
                h = tf[iso][(optics, "hulthen")]
                v = tf[iso][(optics, "vmc")]
                print(f"{iso:<8}{optics:<22}{h[3]:<26}{h[0]:16.4f}{v[0]:12.4f}"
                      f"{v[0]/max(h[0], 1e-12):9.3f}")
        print()
        print("main window / near-beam tail split:")
        for iso in ("6Li", "7Li"):
            for optics in OPTICS:
                for wave in WAVES:
                    a, mn, nb, name = tf[iso][(optics, wave)]
                    print(f"  {iso} {optics:<20} {wave:<8} all {a:.4f} = main "
                          f"{mn:.4f} + near-beam {nb:.4f}")
    # write_markdown() documents the 10 x 99.5 GeV/u point only (config 1),
    # regardless of --configs, so its format stays unchanged bit for bit.
    tf = tf_by_cfg.get(1, tf_by_cfg[configs[0]])

    # ---- A_zz^tag --------------------------------------------------------
    print()
    print("=" * 92)
    print("3.  TAGGED TENSOR ASYMMETRY A_zz^tag(k), 6Li, YR high-acceptance")
    print("=" * 92)
    az = azz_curves(lp)
    order = ("hulthen", "hulthen-P_D_vmc", "vmc", "vmc-flipD")
    kh = az["hulthen"][0]
    print("NaN below k ~ 0.17 GeV: the Roman-Pot envelope accepts nothing "
          "there, so the\nacceptance-weighted numerator AND denominator are "
          "both zero.")
    print(f"{'k [GeV]':>9}" + "".join(n.rjust(18) for n in order)
          + f"{'vmc - hulthen':>15}")
    rows = []
    for k0 in K_ROWS:
        i = int(np.argmin(np.abs(kh - k0)))
        vals = [az[n][1][i] for n in order]
        rows.append((kh[i], vals))
        print(f"{kh[i]:9.4f}" + "".join(f"{v:18.4f}" for v in vals)
              + f"{vals[2] - vals[0]:15.4f}")

    if args.markdown:
        write_markdown(args.markdown, moments, tf, rows, args.events, args.seed)
        print(f"\nAppended to {args.markdown}")


def write_markdown(path, moments, tf, rows, n_events, seed):
    text = path.read_text()
    if MARKER in text:
        text = text[:text.index(MARKER) + len(MARKER)] + "\n"
    else:
        text = text.rstrip() + "\n\n" + MARKER + "\n"

    body = ["", "## Impact on the tagged pipeline",
            "",
            "Generated by `validation/vmc_tag_fractions.py`.  "
            f"{n_events} events per point, e × ion = 10 × 99.5 GeV/u "
            f"(`--config 1`), seed {seed}; `tensor-thirds` plan for ⁶Li, "
            "`helicity-flip` for ⁷Li (J = 3/2 has no equal-thirds pattern).",
            "",
            "### Moments on the `TaggedModel` grid (k ∈ [1e-4, 1.2] GeV, 280 cells)",
            "",
            "| channel | ⟨k⟩ [GeV] | P(k>0.2) | P(k>0.3) | P(k>0.45) | P_D |",
            "|---|---|---|---|---|---|"]
    for (iso, wave), (mk, t2, t3, t45, pd) in moments.items():
        pd_s = "—" if np.isnan(pd) else f"{pd:.5f}"
        body.append(f"| {iso} {wave} | {mk:.4f} | {t2:.4f} | {t3:.4f} "
                    f"| {t45:.4f} | {pd_s} |")

    body += ["", "### Roman-Pot tag fractions", "",
             "| isotope | optics | envelope | Hulthen β=0.30 | VMC AV18 | ratio |",
             "|---|---|---|---|---|---|"]
    for iso in ("6Li", "7Li"):
        for optics in OPTICS:
            h = tf[iso][(optics, "hulthen")]
            v = tf[iso][(optics, "vmc")]
            body.append(f"| {iso} | {optics} | `{h[3]}` | {h[0]:.4f} "
                        f"| {v[0]:.4f} | {v[0]/max(h[0],1e-12):.3f} |")

    body += ["", "### Tagged tensor asymmetry A_zz^tag(k), ⁶Li, YR high-acceptance", "",
             "`hulthen-P_D_vmc` is the analytic shapes at the VMC D-state "
             "probability (it isolates the P_D change from the shape change); "
             "`vmc-flipD` is the VMC channel with the D-wave table negated (the "
             "control on the S–D relative sign, the one input a momentum "
             "density cannot supply).  NaN below k ≈ 0.17 GeV: the Roman-Pot "
             "envelope accepts nothing there.",
             "",
             "| k [GeV] | hulthen β=0.30 | hulthen P_D=VMC | VMC AV18 "
             "| VMC D-sign flipped | VMC − hulthen |",
             "|---|---|---|---|---|---|"]
    for k, vals in rows:
        body.append("| %.4f | " % k
                    + " | ".join("%.4f" % v for v in vals)
                    + " | %+.4f |" % (vals[2] - vals[0]))
    body.append("")
    path.write_text(text + "\n".join(body) + "\n")


if __name__ == "__main__":
    main()
