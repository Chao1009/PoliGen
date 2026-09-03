"""`lipolgen-run` -- one command line for a whole LiPolGen run.

    lipolgen-run --isotope 6Li --config 1 --channel tagged-alpha \
                 --plan tensor-thirds --events 100000 --seed 1 \
                 --hepmc out.hepmc --npz out.npz [--hadronize]

Every switch has a config-file twin: `--config-file run.json` (or run.yaml,
if PyYAML is installed) is read into the same keyword set, and explicit
command-line switches win over the file.  The file is a flat mapping of the
long option names with the dashes turned into underscores:

    {"isotope": "6Li", "config": 1, "channel": "tagged-alpha",
     "plan": "tensor-thirds", "events": 100000, "seed": 1,
     "optics": "yr-high-acceptance", "pzz": 0.6,
     "coherent": {"f0": 0.04, "slope_b": 50.0}}

Outputs
    --npz       the columnar sample (`Pipeline.generate()`), one array per
                column plus a JSON `meta` -- `lipolgen.export.write_columns_npz`
    --hfs-npz   the polligen `HFSSample` .npz of the hadronic final state;
                needs --hadronize (there are no T2 hadrons without it)
    --hepmc     HepMC3 Asciiv3, the full event record
"""

import argparse
import json
import os
import sys
import time

import numpy as np

from . import _lipolgen as _l
from . import export
from . import (CHANNELS, CLUSTER_WAVES, FSI, OPTICS, PLANS, RC, TRITON_SFS,
               ion_spin,
               make_config, make_plan)


def _load_config_file(path):
    with open(path) as f:
        text = f.read()
    if path.lower().endswith((".yaml", ".yml")):
        try:
            import yaml
        except ImportError:
            raise SystemExit(
                "%s looks like YAML but PyYAML is not installed; use JSON "
                "or `pip install pyyaml`" % path)
        data = yaml.safe_load(text)
    else:
        data = json.loads(text)
    if not isinstance(data, dict):
        raise SystemExit("%s must hold a mapping of option names" % path)
    return data


def build_parser():
    p = argparse.ArgumentParser(
        prog="lipolgen-run",
        description="Generate a LiPolGen run and write it out.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--config-file", default=None,
                   help="JSON (or YAML, if PyYAML is installed) file of the "
                        "same options; command-line switches win")
    p.add_argument("--isotope", default=None, help="6Li, 7Li, d")
    p.add_argument("--config", type=int, default=None,
                   help="beam configuration index 0/1/2 = low/mid/top")
    p.add_argument("--channel", default=None, choices=sorted(CHANNELS),
                   help="physics channel")
    p.add_argument("--plan", default=None, choices=sorted(set(PLANS)),
                   help="spin run plan")
    p.add_argument("--events", type=int, default=None,
                   help="fixed event count (exclusive with --lumi)")
    p.add_argument("--lumi", type=float, default=None,
                   help="integrated luminosity [pb^-1] (exclusive with "
                        "--events)")
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--run", type=int, default=None, help="run number")
    p.add_argument("--optics", default=None, choices=sorted(OPTICS),
                   help="far-forward envelope the route label is priced at")
    p.add_argument("--pz", type=float, default=None, help="fill P_z")
    p.add_argument("--pzz", type=float, default=None, help="fill P_zz")
    p.add_argument("--pe", type=float, default=None,
                   help="electron polarization")
    p.add_argument("--rel-lumi-offset", type=float, default=None,
                   help="relative-luminosity offset on the plan's own category")
    p.add_argument("--cluster-beta", type=float, default=None,
                   help="short-range scale of the cluster radial waves")
    p.add_argument("--p-d", type=float, default=None,
                   help="D-state probability (6Li alpha tag)")
    p.add_argument("--cluster-wave", choices=sorted(CLUSTER_WAVES),
                   default=None,
                   help="cluster radial forms: 'hulthen' (default, "
                        "bit-compatible) or 'vmc' (ANL VMC tables; ignores "
                        "--cluster-beta and --p-d on the lithium alpha tags)")
    p.add_argument("--triton-sf", choices=sorted(TRITON_SFS), default=None,
                   help="triton spectral function of the 7Li alpha tag's T1 "
                        "breakup: 'hulthen' (default, the sequential "
                        "two-body decay, bit-compatible) or 'ciofi-simula' "
                        "(the Ciofi degli Atti-Simula n0 + n1 model: "
                        "k-dependent n0/(n0+n1) branching and the third "
                        "channel, struck n -> a (p n) continuum)")
    p.add_argument("--fsi", choices=sorted(FSI), default=None,
                   help="FSI of the DIS debris with the tagged spectator, as "
                        "a per-event WEIGHT on Event.weight (never a shift): "
                        "'off' (default, the plane-wave impulse "
                        "approximation), 'glauber-cluster' (coherent "
                        "X-cluster amplitude, Glauber-shadowed) or "
                        "'glauber-nucleon' (unshadowed A*sigma_XN single-"
                        "scattering bracket).  Tagged channels only")
    p.add_argument("--fsi-sigma-mb", type=float, default=None,
                   help="sigma_XN [mb] the FSI weight is built at; 40 = free "
                        "hadron.  The documented band is 20-40 mb: run BOTH "
                        "ends as a systematic, never pin one row alone")
    p.add_argument("--rc", choices=sorted(RC), default=None,
                   help="tensor-sector radiative corrections as OPT-IN "
                        "WEIGHTS (rc.hpp): 'off' (default) is today bit for "
                        "bit; 'tensor-band' adds rc_tensor_lo/rc_tensor_hi "
                        "(the band on the tensor part of the rate) and "
                        "rc_tail (the 6Li elastic + unpolarised "
                        "quasi-elastic radiative tails).  Never a momentum "
                        "shift; never on Event.weight.  NEEDS AN UNPOLARISED "
                        "BEAM (lam_e*pe = 0) -- the whole A_zz programme "
                        "assumes one -- so on 7Li no --plan this CLI can "
                        "build reaches it and the band is API-only there "
                        "(USAGE sec. 7b).  rc_tail is the t-PEAK ONLY: a "
                        "LOWER BOUND on the dilution, validated for 6Li at "
                        "Q^2 >= 20 GeV^2 and low by ~4x at fixed-target "
                        "kinematics")
    p.add_argument("--rc-delta-low-x", type=float, default=None,
                   help="the low-x band edge (default 0.30, the conservative "
                        "end of Gakh-Shekhovtsova's UNCITED 10-30 %%; 0.19 is "
                        "the residual HERMES actually achieved). BAND IT: "
                        "run both, never quote one row alone")
    p.add_argument("--rc-delta-high-x", type=float, default=None,
                   help="the high-x band edge (default 0.015, E12-13-011)")
    p.add_argument("--rc-fq-scale", type=float, default=None,
                   help="+-100 %% systematic on the 6Li quadrupole form "
                        "factor.  The tensor tail is QUADRATIC in it, so the "
                        "band must be RUN (0.0, 1.0, 2.0), never rescaled "
                        "from one run")
    p.add_argument("--rc-tail-tensor-scale", type=float, default=None,
                   help="multiplier on the 6Li MAGNETIC form factor -- the "
                        "eta*F_m^2 tensor sector, which --rc-fq-scale does "
                        "NOT span (default 1.0; run 0.5 and 2.0)")
    p.add_argument("--rc-qe-suppression", type=float, default=None,
                   help="multiplier on the UNPOLARISED quasi-elastic "
                        "radiative tail, standing in for POLRAD Eq. (44)'s "
                        "S_E/S_M/S_EM factors (default 1.0 = no suppression, "
                        "the conservative direction; run 0.0 and 0.5)")
    p.add_argument("--inclusive-b1", action="store_true", default=None,
                   help="put an inclusive b1 in the struck cluster's kernel")
    p.add_argument("--coherent-f0", type=float, default=None)
    p.add_argument("--coherent-slope-b", type=float, default=None)
    p.add_argument("--coherent-amp", type=float, default=None)
    p.add_argument("--nthreads", type=int, default=None,
                   help="generation threads (forced to 1 with --hadronize)")
    p.add_argument("--hadronize", action="store_true", default=None,
                   help="run the T2 (PYTHIA 8) tier on every event")
    p.add_argument("--coherent-t2", choices=("pomeron", "off"), default=None,
                   help="what --hadronize does with a coherent event: "
                        "'pomeron' (default) hadronizes the gamma*-Pomeron "
                        "system on a PYTHIA Pomeron beam (Beams:idA = 990); "
                        "'off' skips the third PYTHIA instance and leaves "
                        "coherent records at T0 (still conserving)")
    p.add_argument("--pom-set", type=int, default=None,
                   help="PYTHIA PDF:PomSet, the Pomeron parton densities of "
                        "the coherent T2 tier (6 = H1 2006 Fit B LO, "
                        "PYTHIA's own default)")
    p.add_argument("--pom-rescale", type=float, default=None,
                   help="PYTHIA PDF:PomRescale, the overall Pomeron-PDF "
                        "normalization (cancels out of the bridge's "
                        "per-event flavour draw; recorded for "
                        "reproducibility)")
    p.add_argument("--hepmc", default=None, help="HepMC3 Asciiv3 output file")
    p.add_argument("--npz", default=None, help="columnar .npz output file")
    p.add_argument("--hfs-npz", default=None,
                   help="polligen HFSSample .npz (needs --hadronize)")
    p.add_argument("--quiet", action="store_true", default=None)
    return p


#: defaults applied after the config file and the command line are merged
DEFAULTS = dict(isotope="6Li", config=1, channel="inclusive",
                plan="tensor-thirds", events=100000, lumi=0.0, seed=20260713,
                run=1, optics="yr-high-acceptance", pz=0.7, pzz=0.6, pe=0.7,
                rel_lumi_offset=0.0, nthreads=1, hadronize=False,
                quiet=False,
                hepmc=None, npz=None, hfs_npz=None, cluster_beta=None,
                p_d=None, cluster_wave="hulthen", triton_sf="hulthen",
                fsi="off", fsi_sigma_mb=40.0,
                # CONVENTIONS.md "no physics number defined twice": the band
                # anchors live in rc.hpp and are exported, so read them rather
                # than retyping 0.30 / 0.015 here.
                rc="off", rc_delta_low_x=_l.RC_DELTA_LOW_X,
                rc_delta_high_x=_l.RC_DELTA_HIGH_X,
                rc_fq_scale=1.0, rc_tail_tensor_scale=1.0,
                rc_qe_suppression=1.0,
                coherent_t2="pomeron", pom_set=6, pom_rescale=1.0,
                inclusive_b1=False, coherent=None)


def resolve(argv=None):
    """Merge defaults < config file < command line into one option dict."""
    parser = build_parser()
    args = vars(parser.parse_args(argv))
    opts = dict(DEFAULTS)
    if args.get("config_file"):
        fromfile = _load_config_file(args["config_file"])
        unknown = set(fromfile) - set(DEFAULTS) - {"coherent"}
        if unknown:
            raise SystemExit("unknown option(s) in %s: %s"
                             % (args["config_file"], ", ".join(sorted(unknown))))
        opts.update(fromfile)
    coh = dict(opts.get("coherent") or {})
    for k in ("f0", "slope_b", "amp"):
        v = args.pop("coherent_" + k, None)
        if v is not None:
            coh[k] = v
    args.pop("config_file", None)
    for k, v in args.items():
        if v is not None:
            opts[k] = v
    opts["coherent"] = coh or None
    return opts


def main(argv=None):
    opts = resolve(argv)
    quiet = bool(opts["quiet"])
    say = (lambda *a: None) if quiet else (lambda *a: print(*a))

    if opts["events"] and opts["lumi"]:
        raise SystemExit("--events and --lumi are exclusive")
    cfg = make_config(isotope=opts["isotope"], config=opts["config"],
                      channel=opts["channel"], events=opts["events"],
                      lumi_pb=opts["lumi"], seed=opts["seed"], run=opts["run"],
                      optics=opts["optics"], cluster_beta=opts["cluster_beta"],
                      p_d=opts["p_d"], cluster_wave=opts["cluster_wave"],
                      triton_sf=opts["triton_sf"],
                      fsi=opts["fsi"], fsi_sigma_mb=opts["fsi_sigma_mb"],
                      rc=opts["rc"],
                      rc_delta_low_x=opts["rc_delta_low_x"],
                      rc_delta_high_x=opts["rc_delta_high_x"],
                      rc_fq_scale=opts["rc_fq_scale"],
                      rc_tail_tensor_scale=opts["rc_tail_tensor_scale"],
                      rc_qe_suppression=opts["rc_qe_suppression"],
                      inclusive_b1=opts["inclusive_b1"],
                      coherent=opts["coherent"])
    plan = make_plan(opts["plan"], j=ion_spin(cfg.isotope), pz=opts["pz"],
                     pzz=opts["pzz"], pe=opts["pe"],
                     rel_lumi_offset=opts["rel_lumi_offset"])

    nthreads = max(1, int(opts["nthreads"]))
    bridge = None
    if opts["hadronize"]:
        if not _l.HAVE_PYTHIA8:
            raise SystemExit("this build has no PYTHIA 8 tier")
        beams = _l.default_configs(cfg.isotope)[cfg.beam_config]
        t0 = time.time()
        popts = _l.PythiaBridgeOptions()
        popts.coherent_t2 = (_l.CoherentT2.Off if opts["coherent_t2"] == "off"
                             else _l.CoherentT2.Pomeron)
        popts.pom_set = int(opts["pom_set"])
        popts.pom_rescale = float(opts["pom_rescale"])
        bridge = _l.PythiaBridge(beams, popts)
        say("PYTHIA 8 bridge ready in %.1f s" % (time.time() - t0))
        _l.set_pythia_hadronizer(cfg, bridge)
        nthreads = 1                      # PythiaBridge is not re-entrant
    if opts["hfs_npz"] and not opts["hadronize"]:
        raise SystemExit("--hfs-npz needs --hadronize (no T2 hadrons "
                         "otherwise)")

    t0 = time.time()
    p = _l.Pipeline(cfg, plan)
    t_setup = time.time() - t0
    say("%s  %s  %s" % (_l.pipeline_channel_name(cfg.channel),
                        p.beam_config.label(), p.optics.name))
    say("  setup %.2f s, %d events over %d categories, sigma = %.6g pb"
        % (t_setup, p.size(), len(plan), p.sigma_pb()))
    for name, sig, cnt in zip([c.name for c in plan.categories],
                              p.sigma_per_category_pb(), p.counts()):
        say("    %-10s sigma = %12.6g pb   N = %d" % (name, sig, cnt))
    if p.fsi_weight is not None:
        fw = p.fsi_weight
        say("  FSI %s: sigma_XN = %g mb (band 20-40; never quote one row "
            "alone), sigma_X%s = %.1f mb, survival = %.4f"
            % (_l.pipeline_fsi_name(cfg.fsi), fw.sigma_eff_mb(0.0),
               p.tagged_channel.base.spectator,
               fw.sigma_cluster_mb(fw.sigma_eff_mb(0.0)), fw.survival()))
    if p.rc_model is not None:
        r = p.rc_model
        say("  RC %s: delta(x) = %.4g at x=0.01, %.4g at x=0.1; band %s; "
            "tail %s"
            % (_l.rc_mode_name(cfg.rc), r.delta(0.01), r.delta(0.1),
               "on" if r.applies else "OFF",
               "on" if r.tail_applies else "OFF"))
        # Print the exclusion reason ONCE, and only the lines that describe
        # what the run actually did: `ff_provenance` embeds the same reason
        # when no form factor was built, and the tail statistics are
        # meaningless on a channel with no tail.
        if not (r.applies and r.tail_applies):
            say("     %s" % r.exclusion_reason)
        if r.tail_applies:
            say("     6Li FF: %s" % r.ff_provenance)
            say("     clipped %.3g%% of tail NODES (by y-band %.3g / %.3g / "
                "%.3g%%; the y > 0.9 one is where Y+ ~ 1/(1-y) bites)"
                % ((100.0 * r.clipped_cell_fraction,)
                   + tuple(100.0 * f for f in r.clipped_fraction_by_y)))
            say("     tails: elastic + UNPOLARISED quasi-elastic (Eq. 44, "
                "Pauli S(q) at k_F = %.3g GeV x %.2g); the POLARISED QE tail "
                "is NOT priced (Zhou et al., PRL 82 (1999) 687)"
                % (cfg.rc_options.qe_kf_gev, cfg.rc_options.qe_suppression))
            say("     the tail is the t-PEAK ONLY (POLRAD Eqs. 37-39): its "
                "absolute normalisation is NOT validated against an exact "
                "tail.  Measured (test_rc.cpp T8(c)) against the leading-log "
                "s-/p-peaks: > 99 % of the total for 6Li at Q^2 >= 20 "
                "GeV^2, but only 23 % at the HERMES deuteron point.  Read "
                "rc_tail as a LOWER BOUND on the dilution.")
        if r.applies:
            say("     the weights are on Event.rc_weights and NOT on "
                "Event.weight: a systematic variation and a background, not "
                "a correction.")
            say("     Gakh-Shekhovtsova hep-ph/0403262 has ZERO INSPIRE "
                "citations: band it (--rc-delta-low-x 0.19 / 0.30, "
                "--rc-fq-scale 0 / 1 / 2, --rc-qe-suppression 0 / 0.5 / 1), "
                "never quote one row alone.")

    # Events are only materialized when something needs the records: the HFS
    # exporter always, HepMC only when the T2 tier is on (regenerating a
    # hadronized run to stream it out would double the cost; a bare T0 run is
    # cheap enough to regenerate inside `Pipeline.write_hepmc`, which streams
    # and stores nothing).
    need_events = bool(opts["hfs_npz"] or (opts["hepmc"] and opts["hadronize"]))
    t0 = time.time()
    cols = p.generate(0, need_events, nthreads)
    dt = time.time() - t0
    n = int(cols["x"].size)
    say("  generated %d events in %.3f s = %.4g ev/s (%d thread%s)"
        % (n, dt, n / max(dt, 1e-12), nthreads, "" if nthreads == 1 else "s"))
    if p.rc_model is not None and p.rc_model.applies:
        # EVENT-level clipping, which the node fractions above are not: nodes
        # are not event-weighted, and the node statistic is computed at
        # q_n = 0 while the per-event clip carries the (q_n/6) tensor term.
        m = cols["meta"]
        say("     RC clipped %d / %d events on the tail (%.3g%%) and %d "
            "(%.3g%%) on the band (Event.rc_clipped carries the bits)"
            % (m["rc_clipped_tail_events"], n,
               100.0 * m["rc_clipped_tail_event_fraction"],
               m["rc_clipped_band_events"],
               100.0 * m["rc_clipped_band_event_fraction"]))
    if bridge is not None:
        st = bridge.stats
        say("  PYTHIA: %d ok, %d failed, %d retries" %
            (st.n_ok, st.n_failed, st.n_retries))

    if _l.is_tagged(cfg.channel) or cfg.channel == _l.PipelineChannel.CoherentLi6:
        route = cols["route"]
        say("  tag fraction at %s: %.4f (main %.4f, near-beam %.4f)"
            % (p.optics.name, float(np.mean(export.rp_accepted(cols))),
               float(np.mean(route == _l.Route.RomanPots)),
               float(np.mean(route == _l.Route.RPNearBeam))))

    if opts["npz"]:
        export.write_columns_npz(cols, opts["npz"])
        say("  wrote %s (%.1f MB)"
            % (opts["npz"], os.path.getsize(opts["npz"]) / 1e6))
    if opts["hfs_npz"]:
        meta = dict(cols["meta"])
        export.write_hfs_npz(cols["events"], opts["hfs_npz"], meta=meta)
        say("  wrote %s (%.1f MB)"
            % (opts["hfs_npz"], os.path.getsize(opts["hfs_npz"]) / 1e6))
    if opts["hepmc"]:
        if not _l.HAVE_HEPMC3:
            raise SystemExit("this build has no HepMC3 writer")
        t0 = time.time()
        if "events" in cols:
            with _l.HepMC3Writer(opts["hepmc"]) as w:
                for ev in cols["events"]:
                    w.write(ev)
        else:
            p.write_hepmc(opts["hepmc"])
        say("  wrote %s (%.1f MB) in %.1f s"
            % (opts["hepmc"], os.path.getsize(opts["hepmc"]) / 1e6,
               time.time() - t0))
    return 0


if __name__ == "__main__":
    sys.exit(main())
