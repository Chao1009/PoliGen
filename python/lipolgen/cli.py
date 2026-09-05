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
from . import (B1_MODELS, B1_UNPOL, CHANNELS, CLUSTER_WAVES, FSI, OPTICS,
               PLANS, RC, TRITON_SFS,
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
                        "--cluster-beta and --p-d on the lithium alpha tags).  "
                        "TAGGED CHANNELS ONLY -- refused on --channel "
                        "inclusive / coherent, where it is never read (their "
                        "6Li constants do not follow it, 11.61 %% away).  On "
                        "tagged-alpha it also selects the EMBEDDED deuteron "
                        "(AV18 fdeut, P_D = 0.0576), which is -2.03 %% on the "
                        "polarised observables against what that path gave "
                        "before 2026-09-04")
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
                        "assumes one -- so on 7Li the ONLY CLI route to the "
                        "band is '--plan helicity-flip --pe 0'; every other "
                        "CLI plan is refused (helicity plans at pe != 0 by "
                        "RcModel, spin-1 tensor plans by Pipeline), and an "
                        "explicit (P_z, T) J = 3/2 fill still needs the API "
                        "(USAGE sec. 7b).  rc_tail is the t-PEAK ONLY: a "
                        "LOWER BOUND on the dilution.  For 6Li at "
                        "Q^2 >= 20 GeV^2, y <= 0.9 it carries the whole tail "
                        "to +0.61 %% of the EVENT-WEIGHTED mean, but PER "
                        "CELL 24.4 %% of that window is off by > 1 %% and the "
                        "worst cell by x6444 (high x, y ~ 0.009); per-cell "
                        "agreement holds only for 0.15 <= y <= 0.7.  Low by "
                        "~4x at fixed-target kinematics")
    p.add_argument("--rc-delta-low-x", type=float, default=None,
                   help="the low-x band edge (default 0.30, the conservative "
                        "end of Gakh-Shekhovtsova's UNCITED 10-30 %%; 0.19 is "
                        "the residual HERMES actually achieved). BAND IT: "
                        "run both, never quote one row alone")
    p.add_argument("--rc-delta-high-x", type=float, default=None,
                   help="the high-x band edge (default 0.015, E12-13-011)")
    p.add_argument("--rc-a-transfer-frac", type=float, default=None,
                   help="price the A = 2 -> A = 6 TRANSFER of the band "
                        "(design Q8).  EVERY number delta(x) interpolates "
                        "between is a DEUTERON number -- HERMES's measured "
                        "low-x residual, Gakh-Shekhovtsova's 10-30 %%, "
                        "E12-13-011's 1.5 %% -- and no A > 2 tensor RC "
                        "calculation exists at all.  This adds "
                        "frac*delta(x) in QUADRATURE, i.e. widens the band by "
                        "sqrt(1 + frac^2): 0.0 (default) = the shipped band, "
                        "which ASSUMES the deuteron fraction transfers "
                        "exactly; 0.5 = known to 50 %% of itself; 1.0 = as "
                        "uncertain as it is large.  There is no measurement "
                        "to prefer any of them -- this PRICES the assumption, "
                        "it does not correct it.  RUN THE EDGES")
    p.add_argument("--rc-fq-scale", type=float, default=None,
                   help="+-100 %% systematic on the 6Li quadrupole form "
                        "factor.  The tensor tail is QUADRATIC in it, so the "
                        "band must be RUN (0.0, 1.0, 2.0), never rescaled "
                        "from one run")
    p.add_argument("--rc-c0-shape", choices=("ho", "vmc-ft"), default=None,
                   help="which 6Li C0 (monopole) shape F_c AND F_q share: "
                        "\"ho\" (default, the unfitted harmonic oscillator "
                        "every published number was made with) or "
                        "\"vmc-ft\", the j0 transform of the committed ANL "
                        "VMC point-proton density at the SAME measured "
                        "<r^2>_point.  A SHAPE, not a multiplier: RUN BOTH "
                        "EDGES.  The tensor fraction of the elastic tail at "
                        "x = 0.1 is +1.6e-4 on one edge and -4.1e-5 on the "
                        "other, i.e. its sign is a band edge and not a result")
    p.add_argument("--rc-tail-model", choices=("t-peak", "t-peak+ll"),
                   default=None,
                   help="which radiative-tail formulation rc_tail carries: "
                        "\"t-peak\" (default, bit for bit every published "
                        "number -- POLRAD Eqs. (37)-(39), (43), the t-peak "
                        "ALONE and therefore a LOWER BOUND) or "
                        "\"t-peak+ll\", which adds the leading-log s- and "
                        "p-peaks to the UNPOLARISED tail.  A STATED MODEL of "
                        "mixed approximation orders (an eta_A quadrature plus "
                        "a single-z leading log), accurate to ~5-10 %%, with "
                        "NO tensor s/p partner -- so it LOWERS the tensor "
                        "fraction of the tail where it bites.  For 6Li at "
                        "EIC Q^2 >= 20 GeV^2, y <= 0.9 it moves the "
                        "EVENT-WEIGHTED mean tail by only +0.61 %%, but it is "
                        "NOT negligible cell by cell: 331 of that window's "
                        "1356 accepted cells (24.4 %%, 28.2 %% of its cross "
                        "section) move by > 1 %% and the worst by x6444, at "
                        "high x and y ~ 0.009 where the soft 1/(1-z) radiator "
                        "D(z_s) reaches 11.5 and the model itself breaks "
                        "down.  A factor ~4 at fixed-target kinematics.  RUN "
                        "BOTH")
    p.add_argument("--rc-tail-tensor-scale", type=float, default=None,
                   help="multiplier on the 6Li MAGNETIC form factor -- the "
                        "eta*F_m^2 tensor sector, which --rc-fq-scale does "
                        "NOT span (default 1.0; run 0.5 and 2.0)")
    p.add_argument("--rc-qe-suppression", type=float, default=None,
                   help="flat multiplier on the UNPOLARISED quasi-elastic "
                        "radiative tail, applied ON TOP of the de Forest-"
                        "Walecka Pauli suppression S(q), which is already on "
                        "by default at 6Li's measured k_F = RC_QE_KF_GEV "
                        "(RcOptions.qe_kf_gev, API only; 0 disables it).  "
                        "This flag is the band knob, not the physics: "
                        "default 1.0 = no EXTRA suppression, the "
                        "conservative direction; run 0.0 / 0.5 / 1.0")
    p.add_argument("--rc-qe-tensor-scale", type=float, default=None,
                   help="price the POLARISED quasi-elastic tail, which is "
                        "otherwise treated as exactly tensor-blind on the "
                        "piece that is 73 %% of rc_tail at x = 0.1 and 99.9 "
                        "%% at x = 0.3.  Default 0.0 = the shipped, "
                        "bit-for-bit tensor-blind tail.  At 1.0 the "
                        "quasi-elastic tail is given the ELASTIC tail's own "
                        "tensor fraction -- a BORROWED magnitude (a coherent "
                        "nuclear quadrupole fraction on an incoherent nucleon "
                        "process), NOT a derived bound, and possibly ~1e2 too "
                        "SMALL at x <= 0.1 because 6Li's elastic tensor "
                        "fraction is anomalously suppressed for a reason the "
                        "quasi-elastic piece does not share.  LINEAR, so one "
                        "run rescales.  Read the magnitude, never the sign")
    p.add_argument("--inclusive-b1", action="store_true", default=None,
                   help="put an inclusive b1 in the struck cluster's kernel")
    p.add_argument("--b1-model", choices=sorted(B1_MODELS), default=None,
                   help="b1 backend of the INCLUSIVE kernel's rank-2 slot: "
                        "'miller' (default; Li6B1(MillerB1) through "
                        "LI6_B1_RANK2_TRANSFER -- bit-for-bit what every "
                        "published inclusive tensor number was made with, and "
                        "the run PLAN's 'toy'), 'cdks' (the same transfer on "
                        "the other CAMP for b1_d, the digitized CDKS Fig. 4 "
                        "column: |b1| two orders of magnitude smaller below "
                        "x ~ 0.1, comparable or larger above it -- peak "
                        "|x b1| 3.34e-4 against Miller's 4.27e-4 at Q2 = 2.5, "
                        "and 8x LARGER at x = 0.3 with the opposite sign.  "
                        "DOUBLED on 2026-09-03: the CDKS column is already "
                        "per nucleon and is no longer halved) or "
                        "'li6-convolution' (the four-term alpha-d convolution "
                        "of b1_nuclear.hpp).  BOTH opt-in models are "
                        "INCLUSIVE CHANNEL ONLY and 6Li ONLY: on a tagged "
                        "channel the alpha-d density is already in the event "
                        "weight, and elsewhere the flag would reach the "
                        "metadata but not the rate.  li6-convolution "
                        "PASSED its A = 2 magnitude gate on 2026-09-03 -- for "
                        "ONE unpolarised nucleon input, not for the shipped "
                        "default.  Ratio 0.843 against the digitized CDKS "
                        "deuteron peak with CDKS's own MSTW2008 LO at their "
                        "Eq. (21); 0.440 on the DEFAULT toy F2, OUTSIDE the "
                        "gate's [0.5, 2] window.  So quote numbers made with "
                        "--b1-unpol mstw; the default toy backend's are not "
                        "covered by the lift.  With CD-Bonn as well the ratio "
                        "is 1.000338 -- a residual below the error of "
                        "digitizing a published figure, NOT three-digit "
                        "agreement with CDKS.  That gate is A = 2 and says "
                        "nothing about the alpha-d step, so both models stay "
                        "band-only")
    p.add_argument("--b1-band-scale", type=float, default=None,
                   help="the MANDATORY 100 %% band on --b1-model cdks and "
                        "li6-convolution: multiplies the WHOLE b1.  "
                        "Q(6Li) = -0.0818(17) fm^2 (LI6_QUADRUPOLE_FM2, "
                        "TUNL; Pyykko's compilation gives -0.0806(6)) "
                        "against Q_d = +0.2859(3) "
                        "fm^2 -- the alpha-d D wave enters the closest "
                        "measured observable with the OPPOSITE sign to the "
                        "deuteron's and nearly cancels it.  RUN 0 / 1 / 2 and "
                        "quote the envelope; never quote one row alone.  "
                        "REFUSED with --b1-model miller (the band is not "
                        "applied to the published numbers, so recording it "
                        "would claim a variation that did not run)")
    p.add_argument("--b1-unpol", choices=sorted(B1_UNPOL), default=None,
                   help="UNPOLARISED structure-function backend that "
                        "--b1-model li6-convolution folds its own F1 against "
                        "(Li6ConvolutionOptions::unpol): 'toy' (default -- "
                        "the inclusive kernel's own ToyF2, handed over as ONE "
                        "shared object, bit-for-bit what every published "
                        "number was made with), 'mstw' (MSTW2008 LO over "
                        "PYTHIA 8's own pdfdata grid: the PDF CDKS computed "
                        "their b1_d with, and THE CONFIGURATION THE A = 2 "
                        "GATE PASSES ON -- G3b 0.843 against 0.440 on the "
                        "toy.  Needs the optional PYTHIA tier) or 'ct18nlo' "
                        "(the phase-D stand-in, kept selectable so the PDF "
                        "systematic can be quoted; needs the optional LHAPDF "
                        "tier).  MEASURED on Li6ConvolutionB1::b1(x, 2.5): "
                        "mstw/toy = 1.848 / 1.276 / 0.817 at x = 0.10 / 0.30 "
                        "/ 0.50 -- up to a factor 1.85 and not monotone.  It "
                        "moves b1 ONLY: the kernel's own f2_source stays "
                        "ToyF2, so the SPIN-BLIND cell cross section is bit "
                        "for bit and only the tensor shift moves.  REFUSED "
                        "with miller and cdks, "
                        "which never read it, and never silently downgraded "
                        "to the toy in a build without the tier")
    p.add_argument("--b1-alpha-d-dwave-weight", type=float, default=None,
                   help="knob on the alpha-d orbital terms (2d) AND (2a) of "
                        "--b1-model li6-convolution -- they are one physical "
                        "effect and scale together.  The SHAPE variant of the "
                        "band; report 0 / 1 / 2 separately from it.  REFUSED "
                        "with miller and cdks, which never read it")
    p.add_argument("--x-max", type=float, default=None,
                   help="upper x edge of the acceptance window (default 1.0, "
                        "Scenario::x_max).  NEEDED by --b1-model cdks and "
                        "li6-convolution: both carry the CDKS camp's b1_d, a "
                        "Q2 = 2.5 DIGITIZATION with no Q2 evolution, and in "
                        "the topmost default cell (x = 0.955) its b1/F1 "
                        "reaches 6.6 resp. 5.6 -- past the point where "
                        "1 + w_avg stays positive and the sampler refuses the "
                        "run.  Use --x-max 0.95.  Miller's b1 is a ratio "
                        "model and does not need it")
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
                rc_a_transfer_frac=0.0,
                rc_fq_scale=1.0, rc_tail_tensor_scale=1.0,
                rc_qe_suppression=1.0, rc_qe_tensor_scale=0.0,
                rc_c0_shape="ho",
                rc_tail_model="t-peak",
                coherent_t2="pomeron", pom_set=6, pom_rescale=1.0,
                inclusive_b1=False,
                b1_model="miller", b1_band_scale=1.0,
                b1_alpha_d_dwave_weight=1.0, b1_unpol="toy", x_max=None,
                coherent=None)


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


#: What the A = 2 gate's G3b clause measured on each `--b1-unpol` backend, for
#: the run banner.  ONE definition, read from `docs/OPEN_ITEMS_SOLUTIONS.md`
#: sec. 10's clause table; the numbers themselves live there and in
#: `tests/test_b1_nuclear.cpp`, not here.
B1_UNPOL_GATE_ROW = {
    "toy": "its own G3b is 0.440, OUTSIDE the [0.5, 2] window",
    "ct18nlo": "its own G3b is 0.719 -- inside the window, but CT18NLO is the "
               "retired stand-in and not CDKS's own PDF, so it is not the row "
               "the lift was read off",
}


def b1_gate_lines(unpol_name):
    """The A = 2 gate lines of the `--b1-model li6-convolution` run banner.

    A pure function of the unpolarised-backend name, so that both branches
    are testable without building a pipeline.

    WHY THE SECOND LINE EXISTS.  The gate's verdict is a statement about ONE
    unpolarised nucleon input -- MSTW2008 LO at CDKS Eq. (21) -- and the
    SHIPPED DEFAULT is not it.  A banner that prints "PASSED" over a run made
    with a different backend hands the reader a pass their own numbers are not
    covered by, which is the overclaim this line exists to stop.  It says
    which side of the line THIS run is on.
    """
    lines = [
        "A = 2 validation gate: PASSED 2026-09-03 -- for the MSTW2008 LO "
        "nucleon input at CDKS Eq. (21)'s delta-function, where G3b's peak "
        "ratio is 0.843 against the digitized CDKS Fig. 4.  The SHIPPED "
        "DEFAULT toy F2 is 0.440, OUTSIDE the [0.5, 2] window.  With the "
        "CD-Bonn wave function as well it is 1.000338 -- a residual below "
        "the error of digitizing a published figure, NOT three-digit "
        "agreement with CDKS.  It is an A = 2 gate and tests NOTHING about "
        "the alpha-d step -- the band below is not optional.  See "
        "docs/open_items/run_2026-09-03/phase_A_numbers.md.",
    ]
    if unpol_name == "mstw":
        lines.append(
            "^ THIS RUN IS in that configuration (b1 unpol = mstw), so the "
            "lift covers these numbers -- as a band, and with the "
            "configuration quoted.")
    else:
        lines.append(
            "^ THIS RUN IS NOT in that configuration: b1 unpol = %s (%s).  "
            "These numbers are NOT covered by the 2026-09-03 lift; rerun "
            "with --b1-unpol mstw before quoting them as physics."
            % (unpol_name,
               B1_UNPOL_GATE_ROW.get(
                   unpol_name, "not a row the gate was measured on")))
    return lines


def require_b1_unpol_tier(name, have_pythia=None, have_lhapdf=None):
    """Refuse `--b1-unpol` for a backend this build cannot construct.

    The same shape as the `--hadronize` tier check below: say which OPTIONAL
    tier is missing, at the command line, instead of letting the binding's
    RuntimeError come out three frames down.  It is NEVER downgraded to the
    toy -- the toy's A = 2 gate ratio is 0.440 against MSTW's 0.843, so a
    silent fallback would relabel exactly the number the flag exists to fix.
    (A build that HAS the tier but no grid file on disk still fails, one
    frame down, out of `MstwSF` naming the file and the directory.)

    The two flags are arguments so both branches are testable in a build that
    happens to have both tiers.
    """
    have_pythia = _l.HAVE_PYTHIA8 if have_pythia is None else have_pythia
    have_lhapdf = _l.HAVE_LHAPDF if have_lhapdf is None else have_lhapdf
    if name == "mstw" and not have_pythia:
        raise SystemExit(
            "--b1-unpol mstw needs the optional PYTHIA 8 tier: MstwSF reads "
            "PYTHIA's own pdfdata/mstw2008lo.00.dat, and this build has no "
            "PYTHIA tier.  Rebuild with -DLIPOLGEN_WITH_PYTHIA=ON, or run "
            "--b1-unpol toy and say so -- it is not silently substituted, "
            "because the two differ by up to a factor 1.85 on b1.")
    if name == "ct18nlo" and not have_lhapdf:
        raise SystemExit(
            "--b1-unpol ct18nlo needs the optional LHAPDF tier, and this "
            "build has none.  Rebuild with -DLIPOLGEN_WITH_LHAPDF=ON and "
            "install the CT18NLO set, or run --b1-unpol toy and say so.")


def main(argv=None):
    opts = resolve(argv)
    quiet = bool(opts["quiet"])
    say = (lambda *a: None) if quiet else (lambda *a: print(*a))

    if opts["events"] and opts["lumi"]:
        raise SystemExit("--events and --lumi are exclusive")
    require_b1_unpol_tier(opts["b1_unpol"])
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
                      rc_a_transfer_frac=opts["rc_a_transfer_frac"],
                      rc_fq_scale=opts["rc_fq_scale"],
                      rc_tail_tensor_scale=opts["rc_tail_tensor_scale"],
                      rc_qe_suppression=opts["rc_qe_suppression"],
                      rc_qe_tensor_scale=opts["rc_qe_tensor_scale"],
                      rc_c0_shape=opts["rc_c0_shape"],
                      rc_tail_model=opts["rc_tail_model"],
                      inclusive_b1=opts["inclusive_b1"],
                      b1_model=opts["b1_model"],
                      b1_band_scale=opts["b1_band_scale"],
                      b1_alpha_d_dwave_weight=opts["b1_alpha_d_dwave_weight"],
                      b1_unpol=opts["b1_unpol"],
                      coherent=opts["coherent"])
    if opts["x_max"] is not None:
        # `Scenario` comes back BY VALUE from the binding, so it is written on
        # a copy and assigned back -- the `cfg.struck` / `cfg.rc_options`
        # pattern.  Left alone when the switch is absent, so the default run
        # is bit for bit.
        sc = cfg.scenario
        sc.x_max = float(opts["x_max"])
        cfg.scenario = sc
    if (cfg.channel == _l.PipelineChannel.Inclusive
            and cfg.b1_model != _l.B1Model.Miller
            and cfg.scenario.x_max > 0.95):
        # Say WHICH flag fixes it, here, instead of letting the sampler throw
        # "negative phi-averaged density for m=1 at x = 0.955" three frames
        # down -- that message names neither --x-max nor the CDKS table, and
        # `--b1-model cdks` without `--x-max` is otherwise a documented flag
        # whose plain invocation always fails.  See docs/USAGE.md sec. 2a,
        # "The top x cell".
        raise SystemExit(
            "--b1-model %s needs --x-max 0.95: both opt-in backends carry the "
            "CDKS camp's b1_d, a Q2 = 2.5 digitization with no Q2 evolution, "
            "and in the topmost default cell (x = 0.955) b1/F1 reaches 3.3 "
            "(cdks) resp. 5.6 (li6-convolution) -- past where the "
            "phi-averaged density 1 + w_avg stays positive, so InclusiveSampler "
            "refuses the run.  See docs/USAGE.md sec. 2a, 'The top x cell'."
            % _l.b1_model_name(cfg.b1_model))
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
    if (cfg.channel == _l.PipelineChannel.Inclusive
            and cfg.kernel is None
            and cfg.b1_model != _l.B1Model.Miller):
        say("  b1 %s: band scale %g, alpha-d D-wave %g (band 0/1/2; never "
            "quote one row alone)"
            % (_l.b1_model_name(cfg.b1_model), cfg.b1_band_scale,
               cfg.b1_alpha_d_dwave_weight))
        if cfg.b1_model == _l.B1Model.Li6Convolution:
            # WHICH unpolarised PDF the convolution's own F1 came from.
            # Printed on every li6-convolution run, "toy" included: the
            # backends differ by up to a factor 1.85 on b1 and by nothing
            # else in the banner, so a run that does not say which one it
            # used is not reproducible from its own log.
            say("     b1 unpol %s (--b1-unpol; F1 of the convolution only "
                "-- the kernel's own f2_source is ToyF2 on every setting, so "
                "the SPIN-BLIND cell cross section does not move with this "
                "flag; the tensor split does)"
                % _l.b1_unpol_name(cfg.b1_unpol))
            for line in b1_gate_lines(_l.b1_unpol_name(cfg.b1_unpol)):
                say("     " + line)
            say("     4 terms: (1) embedded d S wave, (2d)+(2a) alpha-d "
                "D wave (struck d and struck alpha, one physical effect), "
                "(3) CG depolarization.  N_ad = %.6g suppresses all of them; "
                "the non-alpha-d 18 %% of 6Li is given b1 = 0."
                % _l.VMC_N_ALPHA_D_LI6)
        else:
            say("     Miller (HERMES-like) and CDKS (convolution) are "
                "different CAMPS for b1_d; the library does not adjudicate "
                "between them.  Say which one a plot used.")
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
                "Pauli S(q) at k_F = %.3g GeV x %.2g)"
                % (cfg.rc_options.qe_kf_gev, cfg.rc_options.qe_suppression))
            if cfg.rc_options.qe_tensor_scale == 0.0:
                say("     the POLARISED QE tail is NOT priced (Zhou et al., "
                    "PRL 82 (1999) 687): rc_tail treats as exactly "
                    "tensor-blind the piece that is 73 % of it at x = 0.1 "
                    "and 99.9 % at x = 0.3.  --rc-qe-tensor-scale prices "
                    "that omission with a borrowed magnitude.")
            else:
                say("     the POLARISED QE tail is priced by a STAND-IN at "
                    "qe_tensor_scale = %.4g: it is given the ELASTIC tail's "
                    "own tensor fraction, which is a BORROWED magnitude and "
                    "NOT a derived bound, is possibly ~1e2 too small at "
                    "x <= 0.1, and whose SIGN is meaningless here.  Quote it "
                    "as a price tag on an omission, never as the polarised "
                    "quasi-elastic tail."
                    % cfg.rc_options.qe_tensor_scale)
            if r.options.tail_model == _l.RcTailModel.TPeak:
                say("     the tail is the t-PEAK ONLY (POLRAD Eqs. 37-39): "
                    "its absolute normalisation is NOT validated against an "
                    "exact tail.  Measured (test_rc.cpp T8(c), T8(d)) against "
                    "the leading-log s-/p-peaks, for 6Li at Q^2 >= 20 GeV^2 "
                    "and y <= 0.9: it carries the EVENT-WEIGHTED mean tail to "
                    "0.61 %, but PER CELL only to 0.55 % over 0.15 <= y <= "
                    "0.7 -- 331 of that window's 1356 accepted cells (24.4 %) "
                    "are off by > 1 % and the worst by a factor 6444, at "
                    "x = 0.79, y = 0.0088.  It is 62.8 % of the total at the "
                    "y -> 1 edge (x = 0.01, y = 0.985) and 23 % at the HERMES "
                    "deuteron point.  Read rc_tail as a LOWER BOUND on the "
                    "dilution, and run --rc-tail-model t-peak+ll as the other "
                    "edge.")
            else:
                say("     the tail is t-PEAK + the LEADING-LOG s-/p-PEAKS "
                    "(--rc-tail-model t-peak+ll).  A STATED MODEL of mixed "
                    "approximation orders -- an eta_A quadrature plus a "
                    "single-z leading log, accurate to ~5-10 % and with an "
                    "UNCANCELLED soft 1/(1-z) that overshoots as y -> 0 -- "
                    "NOT a controlled O(alpha) expansion.  It has NO tensor "
                    "s/p partner, so it LOWERS the tensor fraction of the "
                    "tail: the tensor part of those peaks is UNKNOWN, not "
                    "zero.  This is the UPPER edge of a band whose lower "
                    "edge is --rc-tail-model t-peak; quote both.")
        if r.applies:
            say("     the weights are on Event.rc_weights and NOT on "
                "Event.weight: a systematic variation and a background, not "
                "a correction.")
            say("     Gakh-Shekhovtsova hep-ph/0403262 has ZERO INSPIRE "
                "citations: band it (--rc-delta-low-x 0.19 / 0.30, "
                "--rc-fq-scale 0 / 1 / 2, --rc-qe-suppression 0 / 0.5 / 1), "
                "never quote one row alone.")
            if cfg.rc_options.a_transfer_frac == 0.0:
                say("     the A = 2 -> A = 6 TRANSFER of delta(x) is NOT "
                    "priced (design Q8): every band anchor is a DEUTERON "
                    "number and no A > 2 tensor RC calculation exists, so "
                    "this run ASSUMES the deuteron fractional RC transfers "
                    "to 6Li exactly.  --rc-a-transfer-frac prices it.")
            else:
                say("     the A = 2 -> A = 6 transfer of delta(x) is priced "
                    "at a_transfer_frac = %.4g, i.e. the band is widened by "
                    "sqrt(1 + f^2) = %.4g.  A PRICE on design Q8, not a "
                    "correction: there is no A > 2 tensor RC calculation to "
                    "calibrate it against."
                    % (cfg.rc_options.a_transfer_frac,
                       (1.0 + cfg.rc_options.a_transfer_frac ** 2) ** 0.5))

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
