# SPDX-License-Identifier: GPL-3.0-or-later
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
               R_SOURCE, UNPOL_SF, POL_SF,
               PLANS, PZZ_MODES, RC, TRITON_SFS,
               _ION_SPIN, ion_spin,
               make_config, make_plan)

#: b1/F1 in the TOPMOST cell of the default acceptance window, per
#: `--b1-model`.  ONE definition, because it is quoted in four places that
#: drifted apart: this module's `--x-max` help and its `--b1-model` refusal,
#: `python/tests/test_b1_model.py` and `docs/USAGE.md` sec. 2a said 6.6 / 5.6,
#: 3.3 / 5.6, 3.3 / 5.6 and 6.6 / 5.6 respectively, and none of the four was
#: what the code computes.  MEASURED 2026-09-16 off the shipped kernel at the
#: default window's top cell (x = 0.954993, Q2 = 167.34413006, 6Li):
#:
#:     python -c "import lipolgen as lg, numpy as np; _l = lg._lipolgen; \
#:       p = _l.Pipeline(lg.make_config(events=100), \
#:                       lg.tensor_thirds_plan(0.7, 0.6)); \
#:       x = np.asarray(p.dis_sampler.x_cells).max(); \
#:       [print(_l.b1_model_name(m), _l.default_inclusive_kernel( \
#:            _l.ion_by_name('6Li'), m).tables(x, 167.34413006) \
#:            .as_dict()) for m in _l.B1Model.__members__.values()]"
#:
#: The 3.3 was the PRE-2026-09-03 half of the cdks figure and the 6.6 was
#: twice the rounded 3.3 -- rounded arithmetic, not a measurement.
B1_TOP_CELL_X = 0.954993
B1_TOP_CELL_Q2 = 167.34413006
B1_TOP_CELL_B1_OVER_F1 = {"cdks": 6.52, "li6-convolution": 5.91,
                          "miller": 0.145}

#: The sentence those numbers are quoted in, written ONCE.
B1_TOP_CELL_NOTE = (
    "both opt-in backends carry the CDKS camp's b1_d, a Q2 = 2.5 "
    "digitization with no Q2 evolution, and in the topmost default cell "
    "(x = %.3f) b1/F1 reaches %.2f (cdks) resp. %.2f (li6-convolution) -- "
    "past where the phi-averaged density 1 + w_avg stays positive, so "
    "InclusiveSampler refuses the run" % (
        B1_TOP_CELL_X, B1_TOP_CELL_B1_OVER_F1["cdks"],
        B1_TOP_CELL_B1_OVER_F1["li6-convolution"]))


def _load_config_file(path):
    # A MISSING OR MALFORMED --config-file IS A REFUSAL, not a traceback: this
    # function already turned "not a mapping" and "PyYAML missing" into
    # SystemExit and let FileNotFoundError and JSONDecodeError through until
    # 2026-09-16.
    try:
        with open(path) as f:
            text = f.read()
    except OSError as e:
        raise SystemExit("--config-file %s: %s" % (path, e))
    if path.lower().endswith((".yaml", ".yml")):
        try:
            import yaml
        except ImportError:
            raise SystemExit(
                "%s looks like YAML but PyYAML is not installed; use JSON "
                "or `pip install pyyaml`" % path)
        try:
            data = yaml.safe_load(text)
        except yaml.YAMLError as e:
            raise SystemExit("--config-file %s is not valid YAML: %s"
                             % (path, e))
    else:
        try:
            data = json.loads(text)
        except json.JSONDecodeError as e:
            raise SystemExit("--config-file %s is not valid JSON: %s"
                             % (path, e))
    if not isinstance(data, dict):
        raise SystemExit("%s must hold a mapping of option names" % path)
    return data


#: `--isotope` spellings that are not the canonical one -> the key
#: `_ION_SPIN` carries.  Built from `_ION_SPIN` itself, so a species added
#: there is spellable in any case the day it lands and this table cannot
#: drift from that one.
_ISOTOPE_CANON = {name.lower(): name for name in _ION_SPIN}


def build_parser():
    p = argparse.ArgumentParser(
        prog="lipolgen-run",
        description="Generate a LiPolGen run and write it out.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--config-file", default=None,
                   help="JSON (or YAML, if PyYAML is installed) file of the "
                        "same options; command-line switches win")
    # `choices` so an unknown species is argparse's own one-line refusal and
    # not `ion_spin`'s KeyError traceback, and so the help lists what is
    # ACTUALLY accepted: the help read "6Li, 7Li, d" until 2026-09-16 while
    # "p" and "3He" ran.  `type=_isotope_name` runs FIRST and normalises the
    # CASE, because `choices` on its own narrowed the surface instead of
    # describing it -- MEASURED 2026-09-16 against b647c16, where
    # `--isotope 6li --channel tagged-d-p` exits 0 and `--isotope 6li
    # --channel inclusive` exits 1 on `ion_spin`'s KeyError:
    #   * every name `_ION_SPIN` carries, IN ANY CASE (6li, 6LI, D, 3HE, ...),
    #     still runs and still resolves to the isotope it did -- that is what
    #     the normaliser restores, and what test_cli.py pins over
    #     `ISOTOPE_SPELLINGS` (`test_every_isotope_casing_resolves_to_its_canonical_name`, `test_a_lower_case_isotope_still_runs`);
    #   * a spelling `_ION_SPIN` does not carry AT ALL (`zzz`, `6 Li`) now
    #     stops at argparse with exit 2.  It used to exit 1 through the
    #     KeyError refusal on a channel that implies no species, and to exit 0
    #     -- silently, the typed species dropped by `channel_isotope`'s rule
    #     -- on one that does.  THAT run does stop running; it is the point of
    #     the flag, and it is written here rather than claimed away.
    p.add_argument("--isotope", default=None, type=_isotope_name,
                   choices=sorted(_ION_SPIN),
                   help="ion species, in any case (6li = 6Li).  6Li "
                        "(default), 7Li and d are the three the physics "
                        "channels are built for; p and 3He are accepted as "
                        "beam species.  A --channel that implies an isotope "
                        "WINS over this flag (a 6Li alpha tag is a 6Li run "
                        "whatever was typed), and the override is said on the "
                        "banner when a species was actually typed")
    p.add_argument("--config", type=int, default=None,
                   help="beam configuration index 0/1/2 = low/mid/top")
    p.add_argument("--channel", default=None, choices=sorted(CHANNELS),
                   help="physics channel")
    p.add_argument("--plan", default=None, choices=sorted(set(PLANS)),
                   help="spin run plan")
    p.add_argument("--events", type=_nonneg_int, default=None,
                   help="fixed event count (exclusive with --lumi)")
    p.add_argument("--lumi", type=float, default=None,
                   help="integrated luminosity [pb^-1]: it SETS the count "
                        "from the luminosity, so --events is not needed and "
                        "is left at 0 when only this is typed (exclusive "
                        "with an --events given here or in --config-file)")
    p.add_argument("--seed", type=_nonneg_int, default=None)
    p.add_argument("--run", type=_nonneg_int, default=None,
                   help="run number")
    p.add_argument("--optics", default=None, choices=sorted(OPTICS),
                   help="far-forward envelope the route label is priced at")
    p.add_argument("--pz", type=float, default=None, help="fill P_z")
    p.add_argument("--pzz", type=float, default=None,
                   help="fill P_zz (the normalised T at J = 3/2).  READ BY "
                        "--plan tensor-thirds / transverse-tensor / "
                        "tensor-flip ALWAYS, and by helicity-flip ONLY AT "
                        "--pzz-mode typed: at the default --pzz-mode ladder "
                        "that plan builds its fill from the MAX-ENTROPY "
                        "ladder at --pz "
                        "(HelicityFlipOptions::use_explicit_pzz is false) and "
                        "does not read this flag -- the fill's own moments "
                        "are printed in the run banner so the difference is "
                        "never silent")
    p.add_argument("--pzz-mode", choices=sorted(PZZ_MODES), default=None,
                   help="WHICH FILL --plan helicity-flip builds, and read by "
                        "that plan alone.  'ladder' (the DEFAULT) takes the "
                        "max-entropy populations at --pz and does NOT read "
                        "--pzz: it is the branch that was always here, so it "
                        "is bit for bit the pre-flag tree by construction.  "
                        "'typed' honours --pzz through spin1_populations "
                        "(J = 1) / spin32_populations (J = 3/2) and REFUSES a "
                        "value outside the plan's domain with the edge named, "
                        "rather than clamping to it -- at --pz 0.7 the "
                        "J = 3/2 domain is 0.26 <= T <= 0.58, so the CLI's "
                        "own default --pzz 0.6 is outside it by 0.02 and "
                        "--plan helicity-flip --pzz-mode typed is refused "
                        "there.  THE TWO FILLS ARE DIFFERENT PHYSICS, not a "
                        "correction of one another.  MEASURED at the "
                        "standard configuration (inclusive, config 1, "
                        "100000 events, seed 20260713, --pz 0.7 --pe 0.7, "
                        "--pzz 0.5 against the ladder): on 7Li NO observable "
                        "moves -- the rank-2 sector is identically zero and "
                        "both fills honour --pz, so sigma agrees to 1 ulp "
                        "(1.970e-16 relative) and A_par to 6.2e-14 of its own "
                        "statistical error -- while 0.1060 %% of the events "
                        "still land in a different (x, Q2) cell; on 6Li sigma "
                        "moves -0.0023527 %% and A_par -4.27e-6 sigma_stat on "
                        "the shipped --b1-model miller (22-36x less on the "
                        "opt-in b1 models).  What moves on BOTH is the "
                        "RECORDED ALIGNMENT -- T 0.4 -> 0.5 at J = 3/2, "
                        "P_zz 0.409403 -> 0.5 at J = 1 -- which is the "
                        "divisor of every tensor estimator: -20 %% resp. "
                        "-18.1195 %% on delta(A_zz) and delta(cos 2phi) at "
                        "fixed N.  Full tables in "
                        "docs/open_items/run_2026-09-06/phase_A_numbers.md "
                        "sec. A3; the fill's own moments are on the banner's "
                        "`fill` line of every run")
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
                        "before 2026-09-04.  SIGN NOTE (2026-09-06): the "
                        "tagged S-D interference sign was inverted before that "
                        "date -- build_amp2 summed psi_L where the amplitude "
                        "needs phi_L = i^L psi_L -- so any A_zz^tag quoted "
                        "from an older run has the WRONG SIGN.  6Li "
                        "A_zz^tag(k = 0.20 GeV) is now -1.207 (hulthen) and "
                        "-0.519 (vmc), was +0.845 and +0.452.  7Li does not "
                        "move at all; P_D does not move; the two dilutions "
                        "move only in their 6th-7th decimal (<= 1.4e-6 abs; "
                        "0.932495 -> 0.932496 on the d control, a VECTOR one); and "
                        "the SPIN-BLIND rate is unmoved to <= 4e-16 relative "
                        "(1-2 ulp).  A TAG FRACTION IS NOT SPIN-BLIND: for a "
                        "tensor-polarised fill it is the same "
                        "acceptance-weighted integral as A_zz^tag and moves "
                        "with it -- at this CLI\'s own default fill (--pz 0.7 "
                        "max-entropy ladder) the 6Li hulthen model prediction "
                        "goes 0.02703 -> 0.02040, -24.5 %%, and the "
                        "tensor-thirds categories go 0.02813 -> 0.01840 (azz+-) "
                        "and 0.01778 -> 0.03722 (azz0).  Only a spin-blind or "
                        "category-averaged fill is invariant: the uniform-M mix "
                        "fraction is unmoved to 15 digits (0.024675932148828 "
                        "hulthen, 0.033810227625842 vmc) and the equal-thirds "
                        "average reproduces it")
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
                        "RcModel, spin-1 tensor plans by make_plan -- in "
                        "Python, before a Pipeline exists), and an "
                        "explicit (P_z, T) J = 3/2 fill still needs the API "
                        "(USAGE sec. 7b).  At the DEFAULT --rc-tail-model "
                        "t-peak, rc_tail is the t-PEAK ONLY: a LOWER BOUND "
                        "on the dilution -- t-peak+ll and polrad-full change "
                        "exactly that, and the banner says which one ran.  "
                        "For 6Li at "
                        "Q^2 >= 20 GeV^2, y <= 0.9 it carries the whole tail "
                        "to +0.61 %% of the EVENT-WEIGHTED mean at --pz 0 "
                        "(+0.62 %% at this CLI's default P_z = 0.7), but PER "
                        "CELL 24.4 %% of that window is off by > 1 %% and the "
                        "worst cell by x6444 (high x, y ~ 0.009); per-cell "
                        "agreement holds only for 0.15 <= y <= 0.7.  Low by "
                        "~4x at fixed-target kinematics")
    p.add_argument("--rc-delta-low-x", type=float, default=None,
                   help="the low-x band edge (default 0.30, the conservative "
                        "end of Gakh-Shekhovtsova's UNCITED 10-30 %%; 0.19 is "
                        "the residual HERMES actually achieved). BAND IT: "
                        "run both, never quote one row alone. The 0.30 is "
                        "that source's panel value carried UPWARD in x, and "
                        "0.113 is what the panel READS at x = 0.00966, the x "
                        "nearest the 0.01 anchor -- the default is 2.65x it. "
                        "Priced (run 2026-09-06 sec. A2): half-width on A_zz "
                        "at x = 0.01, Q2 = 5 is 1.358e-04 / 1.205e-04 / "
                        "5.117e-05 at 0.30 / 0.266 / 0.113, unmoved at "
                        "x >= 0.16, and the tagged clipped fraction does not "
                        "move at all")
    p.add_argument("--rc-delta-high-x", type=float, default=None,
                   help="the high-x band edge (default 0.015, E12-13-011). "
                        "There is NO alternative to band it against: the "
                        "1.5 %% lives in the UNPUBLISHED proposal only, the "
                        "published companion (arXiv:2506.04506) has no "
                        "radiative-correction discussion at all, and the only "
                        "alternative the record names is an action -- cite it "
                        "by page, or drop the anchor -- not a second value "
                        "(checked 2026-09-06, sec. A2). It pins x >= 0.16 "
                        "exactly, so the low-x anchor decision lives entirely "
                        "below it")
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
    p.add_argument("--rc-tail-model",
                   choices=("t-peak", "t-peak+ll", "polrad-full"),
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
                        "EVENT-WEIGHTED mean tail by only +0.61 %% at --pz 0 "
                        "(+0.62 %% at the default P_z = 0.7), but it is "
                        "NOT negligible cell by cell: 331 of that window's "
                        "1356 accepted cells (24.4 %%, 28.2 %% of its cross "
                        "section) move by > 1 %% and the worst by x6444, at "
                        "high x and y ~ 0.009 where the soft 1/(1-z) radiator "
                        "D(z_s) reaches 11.5 and the model itself breaks "
                        "down.  A factor ~4 at fixed-target kinematics.  RUN "
                        "BOTH.  \"polrad-full\" (2026-09-06) is POLRAD "
                        "Eq. (18) + Appendix B + Eq. (A.4): ONE exact tau_A "
                        "quadrature carrying all THREE peaks WITH the "
                        "s-/p-peaks' own Eq. (A.4) tensor content, so it "
                        "needs no s/p tensor stand-in and --rc-sp-tensor-"
                        "scale is refused on it.  ITS ABSOLUTE NORMALISATION "
                        "IS STILL NOT CHECKED AGAINST MO-TSAI OR ANY "
                        "EXTERNAL EXACT TAIL -- only POLRAD-internally (the "
                        "x_A -> 0 reduction to Eq. (38), unpolarised AND "
                        "tensor) and against the leading-log fallback.  It "
                        "does NOT sit inside the t-peak band: measured over "
                        "the sampler's 3051 accepted cells it lies between "
                        "the two edges on 1725 (56.5 %%), covers a median "
                        "0.6555 of the gap OVER THE 3027 CELLS WHOSE GAP IS "
                        "NONZERO (0.6620 over all 3051; on the other 24 "
                        "t-peak+ll equals t-peak exactly and the fraction is "
                        "undefined), the RATIO of the sigma-weighted mean "
                        "SHIFTS is 0.428 -- a ratio of means, NOT a "
                        "sigma-weighted mean of the per-cell fractions, "
                        "which is 6.483 -- and it "
                        "leaves the band on the far side at the high-x, "
                        "low-y corner -- x4518.3 the t-peak at x = 0.954993, "
                        "Q^2 = 206.68, y = 0.0544 with the ceiling raised.  "
                        "At the SHIPPED rc_options.tail_max = 10 that cell "
                        "and five others clip to x11, which is 1 + tail_max: "
                        "the CEILING, not the model (the 'x11 at x = 0.955, "
                        "Q^2 = 186' this help printed until 2026-09-06 was "
                        "that ceiling mislabelled; that cell is x3449.9 "
                        "unclipped).  Mean rc_tail on 6Li config 1, each "
                        "with its estimator AND its P_z: EVENT-weighted over "
                        "200 000 events at seed 1234 AT THE DEFAULT FILL "
                        "P_z = 0.7, 1.021778527 (t-peak) -> "
                        "1.040274188 (t-peak+ll) -> 1.029702912 "
                        "(polrad-full), and at --pz 0 (the plan both test "
                        "suites use) 1.021836305 -> 1.040368933 -> "
                        "1.029775347 on the same build -- the fill plan "
                        "moves which events the sampler draws, so an "
                        "event-weighted mean is a statement about a "
                        "POLARISATION as well as a seed; "
                        "CELL-cross-section-weighted over the "
                        "3051 accepted cells, 1.02217080 -> 1.04097066 -> "
                        "1.03021634 unclipped, the last becoming 1.03020526 "
                        "at tail_max = 10, and that estimator is P_z-free "
                        "(measured identical on both plans).  "
                        "Costs ~10 s of tail-table build "
                        "and needs rc_options.n_eta >= 64")
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
    p.add_argument("--rc-sp-tensor-scale", type=float, default=None,
                   help="price the TENSOR FRACTION OF THE LEADING-LOG "
                        "s-/p-PEAKS, which --rc-tail-model t-peak+ll "
                        "otherwise adds to the UNPOLARISED numerator alone, "
                        "so that edge LOWERS the tensor fraction of the tail "
                        "by growing the denominator alone (r_T/r_U x0.66 / "
                        "x0.0032 / x6.6e-05 at x = 0.01 / 0.10 / 0.30, "
                        "Q^2 = 5).  NEEDS --rc-tail-model t-peak+ll (refused "
                        "otherwise: the t-peak-only tail computes no s/p "
                        "peaks at all).  Default 0.0.  At 1.0 the ELASTIC "
                        "s-/p-peaks are given the elastic t-peak's own tensor "
                        "fraction.  IT IS A BOUND WITH NO DERIVATION: "
                        "POLRAD's Eq. (38) supplies no tensor s/p peak and "
                        "none is derived here -- but Eq. (18) + Eq. (A.4) "
                        "DOES supply one and --rc-tail-model polrad-full "
                        "computes it, which is why this flag is refused "
                        "there for the OPPOSITE reason (the term ran).  It "
                        "borrows LESS than --rc-qe-tensor-scale "
                        "(one coherent 6Li vertex, two photon topologies) "
                        "but the s/p vertex sits at Q'^2 ~ Q^2 where the "
                        "coherent form factor is DEAD -- so MEASURED the "
                        "shift is EXACTLY 0 at those three points, and "
                        "reaches 583 %% of the band half-width only at "
                        "x = 4.2e-4, y = 0.97 (77 of 3051 cells).  LINEAR, "
                        "so one run rescales.  Read the magnitude, never "
                        "the sign")
    p.add_argument("--inclusive-b1", action="store_true", default=None,
                   help="put an inclusive b1 in the struck cluster's kernel.  "
                        "OFF by default on the tagged channels, deliberately: "
                        "the alpha-d density is already in the event weight "
                        "there, so the per-category TOTALS carry no tensor "
                        "term and the thirds A_zz on sigma_tot is 0 by "
                        "construction (the tagged tensor signal lives in the "
                        "spectator-differential rate).  Turn it on and the "
                        "tagged A_zz moves by the same factors the inclusive "
                        "channel does: -1.557969e-3 (toy) -> -1.951190e-3 "
                        "(ct18nlo, x1.25239) -> -1.963721e-3 (mstw, "
                        "x1.26044).  IT NEEDS --x-max 0.95 OFF THE TOY "
                        "UNPOLARISED BACKEND (measured 2026-09-05, --channel "
                        "tagged-6Li-alpha at the shipped window): the top "
                        "cell x = 0.955, Q2 = 167.3 gives 1 + w_avg = -0.1302 "
                        "under --unpol-sf ct18nlo and -0.03496 under "
                        "--unpol-sf mstw, and InclusiveSampler refuses the "
                        "run by name.  --unpol-sf toy runs unchanged, which "
                        "is why the flag has no unconditional refusal")
    p.add_argument("--b1-model", choices=sorted(B1_MODELS), default=None,
                   help="b1 backend of the INCLUSIVE kernel's rank-2 slot: "
                        "'miller' (default; Li6B1(MillerB1) through "
                        "LI6_B1_RANK2_TRANSFER -- bit-for-bit what every "
                        "published inclusive tensor number was made with, and "
                        "the run PLAN's 'toy'), 'cdks' (the same transfer on "
                        "the other CAMP for b1_d, the digitized CDKS Fig. 4 "
                        "column: |b1| two orders of magnitude smaller below "
                        "x ~ 0.1, comparable or larger above it -- peak "
                        "|x b1| 3.34e-4 at x = 0.766 against Miller's "
                        "4.27e-4 at x = 0.084, both at Q2 = 2.5, and 7.7x "
                        "LARGER at x = 0.3 with the opposite sign.  The two "
                        "peaks sit at DIFFERENT x, so read the comparison off "
                        "the whole digitized range: a scan truncated at "
                        "x = 0.5 finds the CDKS curve's local 5.43e-5 at "
                        "x = 0.332 instead.  DOUBLED on 2026-09-03: the CDKS "
                        "column is already per nucleon and is no longer "
                        "halved) or "
                        "'li6-convolution' (the four-term alpha-d convolution "
                        "of b1_nuclear.hpp).  BOTH opt-in models are "
                        "INCLUSIVE CHANNEL ONLY and 6Li ONLY: on a tagged "
                        "channel the alpha-d density is already in the event "
                        "weight, and validate() REFUSES the flag by name off "
                        "the inclusive channel and off 6Li rather than let it "
                        "reach the metadata but not the rate.  li6-convolution "
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
                        "moves b1 ONLY: the kernel's own f2_source is set "
                        "by a SEPARATE flag, --unpol-sf, and at its shipped "
                        "default 'toy' the SPIN-BLIND cell cross section is "
                        "bit for bit under THIS flag with only the tensor "
                        "shift moving.  The two stay separate because they "
                        "are different quantities -- this one is the "
                        "deuteron F1 inside CDKS Eq. (22), the "
                        "CDKS-comparability choice the A = 2 gate verdict "
                        "rests on -- and validate() refuses the one "
                        "combination that would mislabel them, --b1-unpol "
                        "toy under a non-toy --unpol-sf.  REFUSED "
                        "with miller and cdks, "
                        "which never read it, and never silently downgraded "
                        "to the toy in a build without the tier")
    p.add_argument("--r-source", choices=sorted(R_SOURCE), default=None,
                   help="the ONE R = sigma_L/sigma_T threaded into BOTH "
                        "halves of the 6Li tensor weight.  The shipped "
                        "observable is a RATIO -- A_zz = -(2/3) K/D_phi, "
                        "K = b1 + (1-y)/(x y^2) b2 over "
                        "D_phi = F1 + (1-y)/(x y^2) F2 -- whose numerator's "
                        "R is Li6ConvolutionOptions::r_func and whose "
                        "denominator's is the kernel's own Options::r_func.  "
                        "'unset' (the DEFAULT) installs NEITHER: both hooks "
                        "stay null, resolve_r turns each into r_sigma_lt "
                        "independently, and the run is bit for bit the "
                        "pre-flag tree.  'sigma-lt' installs one r_sigma_lt "
                        "in both and is MEASURED bit-identical to unset "
                        "(sigma_pb, every per-category cross section, every "
                        "generated column, and b1 / K/D_phi / A_zz / the "
                        "cos 2phi amplitude at the six standard points "
                        "x = 0.05/0.10/0.30 x Q2 = 2.5/5, each at y = 0.1, "
                        "0.5 and 0.9) -- "
                        "the wiring's own test, not a physics variation, and "
                        "the provenance table labels it not-read.  'r1998' "
                        "installs the A = 2 gate's SLAC world fit in both.  "
                        "MEASURED at y = 0.5, Q2 = 2.5 against unset: "
                        "K/D_phi (and so A_zz) -3.84 %% / +26.40 %% / "
                        "+3.35 %% and the cos 2phi amplitude -8.33 %% / "
                        "-6.73 %% / -3.14 %% at x = 0.05 / 0.10 / 0.30; the "
                        "weight's shift is y-DEPENDENT (it runs -5.47 %% at "
                        "y = 0.1 to +2.37 %% at y = 0.9 for x = 0.05) "
                        "because only the DENOMINATOR's R moves and D_phi's "
                        "F1 and F2 terms carry different powers of y; the "
                        "numerator-only swap of registry option (ii) is "
                        "y-independent to 2 ulp.  Unlike "
                        "--b1-unpol it moves the UNPOLARISED rate too "
                        "(sigma_pb -0.6847 %% on a 2000-event 6Li run at "
                        "--x-max 0.95), "
                        "because the kernel's r_func reaches F1, F_L, D(y) "
                        "and ToyG1.  REFUSED with miller and cdks, which "
                        "have no numerator R to share")
    p.add_argument("--unpol-sf", choices=sorted(UNPOL_SF), default=None,
                   help="UNPOLARISED structure-function backend of EVERY "
                        "kernel the run builds -- the inclusive kernel, the "
                        "coherent channel that rides its cell cross "
                        "sections (x0.7058 on the shipped 6Li ct18nlo "
                        "coherent run, so it reaches that rate too, which "
                        "--pol-sf does not), and the TAGGED struck-cluster "
                        "kernel, which had no injection point at all before "
                        "2026-09-04.  'toy' (default) is the library's "
                        "ToyF2, labelled TOY in sf.hpp and anchored BY EYE "
                        "('factor-1.5 rate estimates ONLY'); 'mstw' is "
                        "MSTW2008 LO over PYTHIA 8's own pdfdata grid (needs "
                        "the optional PYTHIA tier); 'ct18nlo' is "
                        "LhapdfSF('CT18NLO', 0) (needs the optional LHAPDF "
                        "tier).  MEASURED on 6Li inclusive at config 1: the "
                        "summed accepted cell cross section is 591846.2 pb "
                        "(toy), 472571.9 (ct18nlo, x0.7985), 469556.5 (mstw, "
                        "x0.7934) -- and it is a SHAPE change, F2p at "
                        "Q2 = 10 moving x0.92 / x1.11 / x1.37 / x1.26 "
                        "(ct18nlo) at x = 0.01 / 0.10 / 0.30 / 0.50.  IT IS "
                        "NOT ORTHOGONAL TO --pol-sf: the kernel's default "
                        "ToyG1 is built on the kernel's own UnpolSF, so this "
                        "flag alone moves g1 as well (A1 moves 5-8 %%) -- "
                        "'--pol-sf toy' does NOT mean g1 unchanged.  "
                        "'ct18nlo' puts about a third of the shipped 6Li "
                        "INCLUSIVE window's rate BELOW CT18NLO's own "
                        "Q2 = 1.677 grid floor (36.18 %%; on the COHERENT "
                        "channel the same run's own rate has 44.75 %% below "
                        "it, because f_coh(x) reweights the same cells), "
                        "where LHAPDF extrapolates downward rather than "
                        "freezing; that is not refused, but THIS run's own "
                        "fraction, on THIS run's channel, is printed at the "
                        "banner and recorded in meta.  Never silently "
                        "downgraded to the toy in a build without the "
                        "tier.  With --pol-sf nnpdfpol on '--channel "
                        "tagged-d-p --plan helicity-flip' at --pe != 0 the "
                        "PAIR needs --x-max 0.95 (see that flag)")
    p.add_argument("--pol-sf", choices=sorted(POL_SF), default=None,
                   help="POLARISED structure-function backend of the "
                        "INCLUSIVE and TAGGED kernels -- g1, and through the "
                        "Wandzura-Wilczek relation g2.  IT DOES NOT REACH "
                        "THE COHERENT CHANNEL, and unlike --unpol-sf it "
                        "never could: the coherent rate is spin-independent "
                        "(f_coh(x) times the UNPOLARISED cell cross "
                        "sections) and that channel's tensor signal is the "
                        "recoil azimuth's 1 + c2 cos 2(phi_t - phi_S).  "
                        "MEASURED: every generated column of a coherent run "
                        "is bit-identical between 'toy' and 'nnpdfpol'.  It "
                        "is accepted there rather than refused -- so a "
                        "three-channel scan needs no special case -- and the "
                        "run then RECORDS THE TRUTH: the banner says it did "
                        "not run, meta['pol_sf'] carries 'not read on "
                        "channel coherent-6Li' instead of a backend name, "
                        "and meta['pol_sf_reach'] the reason.  'toy' (default) is "
                        "ToyG1 on the kernel's own UnpolSF; 'nnpdfpol' is "
                        "LhapdfG1('NNPDFpol11_100', 0), the only POLARISED "
                        "set installed here (needs the optional LHAPDF "
                        "tier).  MEASURED on 6Li per-nucleon g1A, which "
                        "A_par tracks to better than 0.1 %% at y = 0.5: "
                        "nnpdfpol/toy = 0.669 / 1.114 / 1.349 / 0.982 / "
                        "0.636 at (x, Q2) = (0.01, 2.5) / (0.10, 10) / "
                        "(0.20, 10) / (0.50, 25) / (0.70, 50) -- x0.64 to "
                        "x1.35 and NOT monotone.  THE NEUTRON IS A SIGN: "
                        "ToyG1's a1n crosses zero near x ~ 0.25 and is "
                        "positive above it while NNPDFpol1.1's g1n stays "
                        "negative to x ~ 0.6, so the shipped toy g1n has the "
                        "WRONG SIGN over roughly 0.25 < x < 0.6 -- which "
                        "matters most on a neutron-tagged run (--channel "
                        "tagged-d-p) and least on isoscalar 6Li.  ON THAT "
                        "VERY CHANNEL 'nnpdfpol' NEEDS --x-max 0.95 when it "
                        "is paired with --unpol-sf ct18nlo|mstw under --plan "
                        "helicity-flip at --pe != 0: the pair's neutron "
                        "g1/F1 reaches 4.40579 in the window's top cell and "
                        "1 + w_avg goes negative, which the sampler refuses "
                        "by name (see --x-max)")
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
                        "the topmost default cell (x = %.3f) its b1/F1 "
                        "reaches %.2f resp. %.2f -- past the point where "
                        "1 + w_avg stays positive and the sampler refuses the "
                        "run.  Use --x-max 0.95.  Miller's b1 is a ratio "
                        "model (%.3f at the same point) and does not need "
                        "it.  ALSO NEEDED, for the "
                        "SAME positivity reason but in the VECTOR sector, by "
                        "'--isotope d --channel tagged-d-p --unpol-sf "
                        "{ct18nlo,mstw} --pol-sf nnpdfpol --plan "
                        "helicity-flip' at any --pe != 0 (the A_par "
                        "measurement on the neutron target, at this CLI's "
                        "default --pe 0.7): the neutron g1/F1 at the same top "
                        "cell (x = 0.955, Q2 = 1119) is 4.40579 on that backend "
                        "PAIR against 0.7228 and 0.3820 on either backend "
                        "alone, giving 1 + w_avg = -0.02414 (ct18nlo) resp. "
                        "-0.1105 (mstw).  The sampler refuses it by name and "
                        "says so; --x-max 0.95 cures it.  AND, third, by "
                        "--inclusive-b1 on a TAGGED channel off the toy "
                        "unpolarised backend: measured 2026-09-05, "
                        "'--channel tagged-6Li-alpha --inclusive-b1 "
                        "--unpol-sf ct18nlo' gives 1 + w_avg = -0.1302 at "
                        "the same top cell (x = 0.955, Q2 = 167.3) and mstw "
                        "gives -0.03496, while --unpol-sf toy runs"
                        % (B1_TOP_CELL_X, B1_TOP_CELL_B1_OVER_F1["cdks"],
                           B1_TOP_CELL_B1_OVER_F1["li6-convolution"],
                           B1_TOP_CELL_B1_OVER_F1["miller"]))
    p.add_argument("--coherent-f0", type=float, default=None)
    p.add_argument("--coherent-slope-b", type=float, default=None)
    p.add_argument("--coherent-amp", type=float, default=None)
    p.add_argument("--coherent-t-max", type=float, default=None,
                   help="|t| ceiling of the coherent channel [GeV^2]; "
                        "default COHERENT_T_MAX_DEFAULT = 0.2.  IT STAYS 0.2 "
                        "because of the ANCHOR RANGE -- the Mantysaari "
                        "deformation table is digitized over |t| <= 0.30 and "
                        "is linear in |t| only as |t| -> 0 -- and NOT "
                        "because of positivity, which binds at 0.245 "
                        "(P_zz = -2) only at the shipped eps_b0 = -0.08 and "
                        "would sit at 2.7990 GeV^2 at the measured 6Li "
                        "quadrupole (derived eps_b0 = -0.0070024; the docs' "
                        "tables say 2.8000, which is arithmetic on the "
                        "ROUNDED -0.0070 -- quote 2.80).  The run banner "
                        "DERIVES it through "
                        "CoherentScenario.t_positivity_edge.  "
                        "Raising it past the positivity edge THROWS.  It "
                        "costs no rate either way: sample_t renormalizes on "
                        "[0, t_max] and exp(-B t_max) = 4.5e-5 of the "
                        "untruncated exponential lies beyond 0.2.  Recorded "
                        "in the npz meta")
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
                        "PYTHIA's own default).  THE COHERENT T2 TIER'S "
                        "LARGEST MODEL SYSTEMATIC, and MEASURED: 20000 "
                        "events per set at 6Li config 1, seed 4242.  It "
                        "moves NOTHING at T0 -- t, x_pom, q2, x and the "
                        "event weight are bit-identical across the 14 sets "
                        "the constructor still admits (1-10, 12-15; set 11 "
                        "is refused, so 14 is the count that reproduces), "
                        "so there is no PomSet band on |t|, x_P, M_X or the "
                        "cross section.  The band is on the HADRONIC FINAL "
                        "STATE, over the 12 genuine DPDF fits (3-10, 12-15) "
                        "about set 6: <n_charged> 3.964 with -2.6%% / +9.9%% "
                        "(min set 9, max set 5), <n_had> 8.521 -2.6%% / "
                        "+10.2%%, <p_T> 0.342 GeV -5.9%% / +10.9%%, and the "
                        "KAON FRACTION 0.0646 with -44%% / +55%%, a factor "
                        "2.8 end to end.  Quote it as an ENVELOPE OVER "
                        "RE-RUNS (one npz per set; no per-event weight maps "
                        "one set onto another).  Sets 1 (toy) and 2 (pi0) "
                        "are not Pomeron fits and are outside the band; set "
                        "11 (PomHISASD) is REFUSED -- it needs setXPom, "
                        "which the bridge never calls, and 100.00%% of its "
                        "events take the e_q^2 fallback.  Set 6 is the only "
                        "LO H1 set, so the band mixes LO and NLO DPDFs: fine "
                        "for an envelope, not for a central value")
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
                pzz_mode="ladder",
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
                rc_sp_tensor_scale=0.0,
                rc_c0_shape="ho",
                rc_tail_model="t-peak",
                coherent_t2="pomeron", pom_set=6, pom_rescale=1.0,
                inclusive_b1=False,
                b1_model="miller", b1_band_scale=1.0,
                b1_alpha_d_dwave_weight=1.0, b1_unpol="toy",
                r_source="unset",
                unpol_sf="toy", pol_sf="toy", x_max=None,
                coherent_t_max=None,
                coherent=None)


def resolve(argv=None):
    """Merge defaults < config file < command line into one option dict."""
    parser = build_parser()
    args = vars(parser.parse_args(argv))
    opts = dict(DEFAULTS)
    fromfile = {}
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
    # LUMINOSITY MODE IS REACHABLE BY TYPING ONE FLAG (2026-09-05).  `DEFAULTS`
    # carries events = 100000, so `--lumi X` alone used to arrive at `main()`
    # with BOTH set and was refused with "--events and --lumi are exclusive" --
    # a mode that only `--events 0 --lumi X` could reach, which no help text
    # and no document said.  Typing `--lumi` and no `--events` anywhere means
    # "take the count from the luminosity"; an `--events` given on the command
    # line or in the config file still collides, and is still refused there.
    if (args.get("lumi") is not None and args.get("events") is None
            and "events" not in fromfile):
        opts["events"] = 0
    # WHO ASKED FOR THE SPECIES.  `DEFAULTS` fills 6Li, so `opts["isotope"]`
    # alone cannot tell a typed species from the filled one -- and `main`'s
    # override notice must fire only for a species the CALLER asked for
    # (2026-09-16: `lipolgen-run --channel tagged-d-p` printed "--isotope 6Li
    # overridden" with no --isotope anywhere on the command line).  The fill
    # happens above and this records what it covered; the command line and
    # the config file both count as typed, the config file taking the same
    # option names the switches do.
    opts["isotope_typed"] = (args.get("isotope") is not None
                             or "isotope" in fromfile)
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


def rank2_zero_lines(report, isotope="7Li", width=74):
    """THE LOUD ZERO: the `--isotope 7Li --channel inclusive` banner block.

    A pure function of the run's own `rank2_input_report` sentence, so it is
    testable without building a pipeline -- `b1_gate_lines`' arrangement.

    THE CLAIM IS NOT WRITTEN HERE.  The first paragraph is the C++
    `rank2_input_report` string verbatim, wrapped; it is the SAME string
    `meta["rank2_input"]` carries, so the banner and the file cannot drift
    apart (docs/CONVENTIONS.md: no physics statement defined twice).  What
    this function adds is the ADVICE, which is a property of the run surface
    and not of the physics.

    Why the advice is what it is, and each clause is measured:

    * The VECTOR sector of a 7Li inclusive run is fine.  `SPIN32_FINITE_GAMMA`
      sec. 4's rank-<=2 theorem plus the measurement in
      `phase_D_li7_rank2.md` sec. 1.6: with a rank-2 slot filled by hand the
      sampler, the tensor weight, the per-category cross sections and the P8
      spin guard are all already correct at J = 3/2.  It is an input that is
      missing, not machinery.
    * The TAGGED 7Li alpha channel is NOT affected: its alignment is in the
      event weight, and it is gated today -- <P2(cos theta_k)> = -T/5 to
      3e-4, `tests/test_tagged.cpp`.
    * A tensor programme on the INCLUSIVE channel needs `--isotope 6Li`,
      where `--b1-model` reaches a real backend.
    """
    import textwrap
    paragraphs = [
        report,
        "^ The rest of this run is unaffected: the unpolarised rate and the "
        "VECTOR sector (g1, A_parallel) are correct at J = 3/2 -- what is "
        "missing is a physics INPUT, not machinery.  Measured with a rank-2 "
        "slot supplied by hand (b1_32 = +0.05*F1, delta_32 = -1e-2*F1, 40000 "
        "events, seed 11): sigma = 565389 / 616516 pb, A_T = -0.043258.",
        "^ For an INCLUSIVE tensor programme use --isotope 6Li, where "
        "--b1-model reaches a real backend.  A TAGGED 7Li run is NOT "
        "affected: --channel tagged-7Li-alpha carries the alpha-t alignment "
        "in the event weight, gated today at <P2(cos theta_k)> = -T/5.",
    ]
    lines = []
    for para in paragraphs:
        for chunk in textwrap.wrap(para, width=width):
            lines.append(("  %s rank-2: " % isotope if not lines else "     ")
                         + chunk)
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


def require_unpol_sf_tier(name, have_pythia=None, have_lhapdf=None):
    """Refuse `--unpol-sf` for a backend this build cannot construct.

    The `require_b1_unpol_tier` shape verbatim: say which OPTIONAL tier is
    missing, at the command line, instead of letting the binding's
    RuntimeError come out three frames down.  It is NEVER downgraded to the
    toy -- the toy is x0.80 on the shipped 6Li rate and up to 24 % away on
    the F2n/F2p the species draws use, so a silent fallback would relabel
    exactly the number the flag exists to fix.
    """
    have_pythia = _l.HAVE_PYTHIA8 if have_pythia is None else have_pythia
    have_lhapdf = _l.HAVE_LHAPDF if have_lhapdf is None else have_lhapdf
    if name == "mstw" and not have_pythia:
        raise SystemExit(
            "--unpol-sf mstw needs the optional PYTHIA 8 tier: MstwSF reads "
            "PYTHIA's own pdfdata/mstw2008lo.00.dat, and this build has no "
            "PYTHIA tier.  Rebuild with -DLIPOLGEN_WITH_PYTHIA=ON, or run "
            "--unpol-sf toy and say so -- it is not silently substituted, "
            "because the two differ by a factor 0.79 on the 6Li rate.")
    if name == "ct18nlo" and not have_lhapdf:
        raise SystemExit(
            "--unpol-sf ct18nlo needs the optional LHAPDF tier, and this "
            "build has none.  Rebuild with -DLIPOLGEN_WITH_LHAPDF=ON and "
            "install the CT18NLO set, or run --unpol-sf toy and say so.")


def require_pol_sf_tier(name, have_lhapdf=None):
    """Refuse `--pol-sf` for a backend this build cannot construct.

    Same rule as `require_unpol_sf_tier`, and the reason is sharper: the
    shipped `ToyG1`'s g1n has the WRONG SIGN over roughly 0.25 < x < 0.6, so
    a silent downgrade would not be a factor, it would be a sign.
    """
    have_lhapdf = _l.HAVE_LHAPDF if have_lhapdf is None else have_lhapdf
    if name == "nnpdfpol" and not have_lhapdf:
        raise SystemExit(
            "--pol-sf nnpdfpol needs the optional LHAPDF tier, and this "
            "build has none.  Rebuild with -DLIPOLGEN_WITH_LHAPDF=ON and "
            "install the NNPDFpol11_100 set (the only polarised set this "
            "tree ships), or run --pol-sf toy and say so -- it is not "
            "silently substituted, because ToyG1's g1n has the wrong SIGN "
            "over roughly 0.25 < x < 0.6.")


def sf_banner_lines(unpol_name, pol_name, q2_min=None, below_frac=None,
                    pol_read=True):
    """The `--unpol-sf` / `--pol-sf` block of the run banner.

    A pure function of the two names, the run's own grid report and the run's
    own polarised-reach report, so every branch is testable without building
    a pipeline -- `b1_gate_lines`'s arrangement, and for the same reason: a
    run that does not say which backend it used is not reproducible from its
    own log.  Printed on EVERY run, the all-default one included, because
    "toy" is a physics choice too and it is the one that costs the most.

    `pol_read` is `_l.pol_sf_is_read(cfg, plan)` -- ONE definition, shared
    with `meta["pol_sf"]` and with the `pol_sf` row of
    `Pipeline.knob_provenance`.  It is read here for the two PRICE clauses
    below, which are true only where g1 actually ran.

    There was a second parameter, `pol_reach` (`_l.pol_sf_reach_report`),
    until 2026-09-16: declared, documented, passed -- and never read in this
    body, so the call that built it was work thrown away on every run.  The
    reach sentence it would have carried is the `pol_sf` row of the KNOB
    PROVENANCE block, which is where reach now lives.

    WHAT THIS BLOCK NO LONGER SAYS, AND WHY.  Until 2026-09-05 its first line
    read "--unpol-sf reaches EVERY kernel this run builds ...; --pol-sf
    reaches the inclusive and the tagged kernels", and a second block said
    "--pol-sf DID NOT RUN HERE" on the one channel that had been looked at.
    Both were hand-written reach claims about two of the run's forty-odd
    knobs, and the first of them was FALSE on the CLI's own default plan (the
    tensor plans carry lam_e = 0, so --pol-sf reaches no kernel there).  The
    reach of every knob is now printed from ONE table --
    `knob_provenance_lines`, off `Pipeline.knob_provenance` -- which cannot
    be true of one knob and stale about the next.  What stays here is what
    the table is not: the measured PRICE of the two backends.
    """
    lines = [
        "  SF backends: unpol %s (--unpol-sf), pol %s (--pol-sf)%s"
        % (unpol_name, pol_name,
           "" if pol_read else "  [--pol-sf did not run: see KNOB PROVENANCE]"),
    ]
    if unpol_name == "toy" and pol_name == "toy":
        lines.append(
            "     BOTH ARE THE TOY BACKENDS.  ToyF2 and ToyG1 are labelled "
            "TOY in sf.hpp and anchored by eye ('factor-1.5 rate estimates "
            "ONLY').  Measured price on 6Li inclusive at config 1: the "
            "accepted cell cross section is x0.80 on ct18nlo and x0.79 on "
            "mstw" + (
                "." if not pol_read else
                ", and ToyG1's g1n has the WRONG SIGN over roughly "
                "0.25 < x < 0.6."))
    if unpol_name != "toy" and pol_name == "toy" and pol_read:
        lines.append(
            "     NOTE: --pol-sf toy does NOT mean g1 is unchanged.  The "
            "kernel's default ToyG1 is built on the kernel's OWN UnpolSF, so "
            "--unpol-sf alone moves g1 too (A1 = g1/F1 moves 5-8 %).  State "
            "both settings next to any number.")
    if unpol_name != "toy":
        # ToyF2 is closed form and has no grid, so this clause is about the
        # selected fit only.  NaN means the backend reports no floor
        # (PYTHIA's MSTWpdf keeps `qsqmin` private; a Custom object reports
        # nothing) -- it is never printed as a 0 that was not measured.
        if below_frac is not None and below_frac == below_frac:   # not NaN
            lines.append(
                "     grid floor Q2 = %.6g; %.2f %% of THIS run's own "
                "accepted RATE is BELOW it, where LHAPDF keeps evolving "
                "downward instead of freezing.  The denominator is the "
                "CHANNEL's own per-cell rate, which on the coherent channel "
                "is sigma_cell x f_coh(x) and not the inclusive cells "
                "(44.75 %% against 36.18 %% on the shipped 6Li ct18nlo run). "
                " There is NO --q2-min flag: raise Scenario.q2_min from "
                "Python to move the window off the extrapolation (it costs "
                "rate -- 472572 -> 301600 pb on the shipped 6Li ct18nlo "
                "inclusive run at q2_min 0.7 -> 1.7)."
                % (q2_min, 100.0 * below_frac))
        else:
            lines.append(
                "     this backend reports no grid floor, so the below-grid "
                "fraction is NOT known and is recorded as NaN, never as 0.")
    return lines


def knob_provenance_lines(rows, width=74):
    """THE run banner's one block for what this run read and what it did not.

    A pure function of `Pipeline.knob_provenance(context)` -- a list of rows
    carrying {name, flag, value, status, reason, label, at_default} -- so it
    is testable without building a pipeline, the `b1_gate_lines` arrangement.

    WHY IT REPLACED THE PER-KNOB LINES.  The banner used to make reach claims
    knob by knob, in prose, wherever somebody had thought to add one: the
    `--unpol-sf` / `--pol-sf` block claimed a reach that was false under the
    default plan, the PomSet block was printed on the coherent channel only
    (so off it the banner was SILENT rather than wrong), and `--cluster-wave`,
    `--triton-sf` and `--inclusive-b1` had no line anywhere.  Five rounds of
    "a knob that did not run was recorded as if it had", each on a new axis.
    One table, printed once, cannot go stale on the axis nobody looked at.

    WHAT IT PRINTS, AND WHAT IT LEAVES TO THE FILE.  Every NOT-READ knob that
    is set away from its default gets its own line WITH THE REASON -- those
    are exactly the ones that mislead, and they are what the five rounds were
    about.  Everything else is named but not explained: the read knobs (values
    for the ones off their defaults), the not-read defaults, and the refused
    axes.  The whole table, reasons included, is in
    `meta["knob_provenance"]`.
    """
    import textwrap
    read = [r for r in rows if r.status == _l.KnobStatus.Read]
    unread = [r for r in rows if r.status == _l.KnobStatus.NotRead]
    refused = [r for r in rows if r.status == _l.KnobStatus.Refused]
    loud = [r for r in unread if not r.at_default]
    quiet = [r for r in unread if r.at_default]
    lines = textwrap.wrap(
        "KNOB PROVENANCE -- what this run READ and what it did not "
        "(%d knobs; meta[\"knob_provenance\"] carries the whole table with "
        "every reason)" % len(rows),
        width=width + 2, initial_indent="  ", subsequent_indent="     ")

    def wrap(head, body):
        return textwrap.wrap(head + body, width=width,
                             initial_indent="     ",
                             subsequent_indent="       ")

    moved = [r for r in read if not r.at_default]
    lines += wrap("read (%d), of which set away from the default: " % len(read),
                  ", ".join("%s = %s" % (r.name, r.value) for r in moved)
                  or "none -- every knob this run reads is at its default")
    if loud:
        lines += textwrap.wrap(
            "NOT READ, and set away from the default -- a value here would "
            "mislead, so the file records the LABEL and the reason:",
            width=width, initial_indent="     ", subsequent_indent="       ")
        for r in loud:
            lines.append("       %s = %s%s"
                         % (r.name, r.value,
                            "  (%s)" % r.flag if r.flag else ""))
            for chunk in textwrap.wrap(r.reason, width=width - 4):
                lines.append("         " + chunk)
    if quiet:
        # GROUPED BY LABEL, not merely counted.  A knob that is not read is
        # still a knob whose value the file records, and the reader has to be
        # able to see WHICH channel, plan or mode put it out of reach -- the
        # `--pol-sf toy` on the coherent channel case, where nothing is
        # promised and "toy" would still name a backend that did not run.
        lines.append("     not read, at their defaults (%d):" % len(quiet))
        groups = {}
        for r in quiet:
            groups.setdefault(r.label or "not read", []).append(r.name)
        for label in sorted(groups):
            lines += textwrap.wrap(
                "%s: %s" % (label, ", ".join(groups[label])), width=width,
                initial_indent="       ", subsequent_indent="         ")
    if refused:
        lines += wrap("refused axes (%d), where validate() throws on any "
                      "other value: " % len(refused),
                      ", ".join(r.name for r in refused))
    return lines


def _isotope_name(text):
    """`--isotope`'s spelling, normalised BEFORE argparse checks `choices`.

    CASE IS THE WHOLE OF IT.  `_ION_SPIN`'s keys are the canonical spellings
    (`3He`, `6Li`, `7Li`, `d`, `p`) and this maps any casing of one of them
    onto its key; anything else is handed back unchanged, so argparse's own
    `choices` refusal names it.  MEASURED 2026-09-16 -- b647c16's
    `python/lipolgen/` on this build's extension module, the two halves of the
    question being argparse's alone: `--isotope 6li --channel tagged-d-p`
    exits 0 there and exited 2 with the bare `choices=` this replaces, while
    `--isotope 6li --channel inclusive` exits 1 there on `ion_spin`'s KeyError,
    exited 2 with the bare `choices=`, exits 0 with this (measured).  A `--channel` that implies its own species drops the typed
    one (`channel_isotope`'s rule), so on that channel the typed spelling
    never reached `ion_spin` and any casing ran.
    """
    return _ISOTOPE_CANON.get(text.lower(), text)


def _nonneg_int(text):
    """A count that reaches an UNSIGNED pybind slot.

    `PipelineConfig::events` / `seed` / `run` are `std::uint64_t`, so a
    negative number met the user as pybind's
    `TypeError: incompatible function arguments ... SupportsInt` with the
    whole overload table pasted under it.  argparse says it in one line
    instead, before anything is built (2026-09-16).
    """
    try:
        v = int(text)
    except ValueError:
        raise argparse.ArgumentTypeError("%r is not an integer" % text)
    if v < 0:
        raise argparse.ArgumentTypeError(
            "%d is negative; this count is unsigned" % v)
    return v


def _refusal(exc):
    """The one-line text of a configuration-time refusal.

    `str(KeyError("6li"))` is `"'6li'"` -- the repr, quotes and all -- which
    is why the message is built from the key rather than from `str`.
    """
    if isinstance(exc, KeyError):
        return ("unknown isotope %s; know %s"
                % (exc.args[0] if exc.args else "?",
                   ", ".join(sorted(_ION_SPIN))))
    return str(exc)


def main(argv=None):
    opts = resolve(argv)
    quiet = bool(opts["quiet"])
    say = (lambda *a: None) if quiet else (lambda *a: print(*a))

    if opts["events"] and opts["lumi"]:
        raise SystemExit("--events and --lumi are exclusive")
    if opts["hepmc"] and not _l.HAVE_HEPMC3:
        # BEFORE the run, not after it.  This check sat below the generation
        # loop until 2026-09-16, so `--events 1000000 --hepmc out.hepmc` on a
        # build without HepMC3 generated the whole run and then refused --
        # while the PYTHIA-tier check a few lines down, and the
        # --hfs-npz/--hadronize check beside it, have always run first.
        raise SystemExit("this build has no HepMC3 writer")
    require_b1_unpol_tier(opts["b1_unpol"])
    require_unpol_sf_tier(opts["unpol_sf"])
    require_pol_sf_tier(opts["pol_sf"])
    # EVERY CONFIGURATION-TIME REFUSAL IS A ONE-LINE MESSAGE, NOT A TRACEBACK
    # -- the rule this module already states for `make_plan` below and for the
    # PYTHIA bridge, and which `make_config` was outside of until 2026-09-16.
    # `PipelineConfig::validate()` throws `RuntimeError` (--fsi
    # glauber-cluster off a tagged channel, --cluster-wave vmc off its
    # channels, --b1-band-scale off {0,1,2}, --rc-sp-tensor-scale on a
    # non-tensor-band --rc, ...), `make_config` itself raises `ValueError`
    # (--config 5), a pybind setter raises `TypeError` (--events -5, --seed -1
    # through the unsigned cast) and `ion_spin` raises `KeyError` (--isotope
    # xx); all four reached the user as a stack trace naming a line of
    # `__init__.py`.  Caught here by TYPE, not by message, so a refusal added
    # to `validate()` later is covered the day it lands.
    try:
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
                          rc_sp_tensor_scale=opts["rc_sp_tensor_scale"],
                          rc_c0_shape=opts["rc_c0_shape"],
                          rc_tail_model=opts["rc_tail_model"],
                          inclusive_b1=opts["inclusive_b1"],
                          b1_model=opts["b1_model"],
                          b1_band_scale=opts["b1_band_scale"],
                          b1_alpha_d_dwave_weight=opts["b1_alpha_d_dwave_weight"],
                          b1_unpol=opts["b1_unpol"],
                          r_source=opts["r_source"],
                          unpol_sf=opts["unpol_sf"], pol_sf=opts["pol_sf"],
                          coherent_t_max=opts["coherent_t_max"],
                          coherent=opts["coherent"])
    except (RuntimeError, ValueError, TypeError, KeyError) as e:
        raise SystemExit(_refusal(e))
    if opts["isotope_typed"] and opts["isotope"] != cfg.isotope:
        # SAY THE OVERRIDE -- AND ONLY WHEN THERE IS ONE TO SAY.
        # `make_config` resolves the isotope a --channel IMPLIES over the one
        # that was typed (`channel_isotope`'s rule), so `--isotope 7Li
        # --channel tagged-d-p` is a deuteron run; until 2026-09-16 nothing on
        # the banner or in the meta said the typed species had been dropped.
        # `opts["isotope_typed"]`, not `opts["isotope"] is not None`: DEFAULTS
        # fills 6Li, so the second test was true on every run and the notice
        # fired for a species nobody typed (`lipolgen-run --channel
        # tagged-d-p` announced "--isotope 6Li overridden").  Typed on the
        # command line or in the config file; filled by DEFAULTS is not typed.
        say("  --isotope %s overridden: --channel %s implies %s"
            % (opts["isotope"], opts["channel"], cfg.isotope))
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
            "--b1-model %s needs --x-max 0.95: %s.  See docs/USAGE.md sec. 2a, "
            "'The top x cell'."
            % (_l.b1_model_name(cfg.b1_model), B1_TOP_CELL_NOTE))
    try:
        plan = make_plan(opts["plan"], j=ion_spin(cfg.isotope), pz=opts["pz"],
                         pzz=opts["pzz"], pe=opts["pe"],
                         rel_lumi_offset=opts["rel_lumi_offset"],
                         pzz_mode=opts["pzz_mode"])
    except (RuntimeError, ValueError, KeyError) as e:
        # A FILL OUTSIDE THE PLAN'S DOMAIN IS A REFUSAL, not a traceback, and
        # it is never CLAMPED to the edge: clamping would publish an alignment
        # nobody typed.  `spin1_populations` / `spin32_populations` name the
        # offending m, its negative population and both edges at this --pz,
        # and `--pzz-mode typed` is what makes that reachable from the
        # command line (--pzz 0.6 at --pz 0.7 is outside the J = 3/2 domain by
        # 0.02: the edge is T = 0.58, where p(-1/2) = 0).  The same clause
        # carries make_plan's spin-1-pattern refusal, which reached the user
        # as a traceback until 2026-09-06, and `ion_spin`'s KeyError on an
        # unknown --isotope until 2026-09-16 (KeyError is in the clause, and
        # `_refusal` strips the quotes `str(KeyError)` would otherwise leave).
        raise SystemExit(_refusal(e))

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
        # ONE unpolarised backend for the whole run (docs/CONVENTIONS.md).
        # The T2 struck-nucleon SPECIES draw is P(p) = Z F2p/(Z F2p + N F2n)
        # on `PythiaBridgeOptions::f2_source`, and the `Pipeline` cannot fill
        # it -- the CLI builds the bridge itself, so `Pipeline` never sees it
        # the way it sees `BreakupOptions::f2`.  Left unset the bridge would
        # keep drawing from ToyF2 while T0 and T1 moved to the selected
        # backend, and the toy's F2n/F2p is up to 24 % away from CT18NLO's at
        # x = 0.5.  None (the toy default) leaves the bridge's own ToyF2 in
        # place, so this line is bit for bit at `--unpol-sf toy`.
        popts.f2_source = cfg.unpol_sf_obj
        try:
            bridge = _l.PythiaBridge(beams, popts)
        except RuntimeError as e:
            # The bridge's own refusals (an unusable `--pom-set`, a beam energy
            # below the nucleon mass, a PYTHIA init failure) are configuration
            # errors typed at the command line, not internal faults: show the
            # message the bridge wrote, not a traceback.  The text is NOT
            # duplicated here -- CONVENTIONS.md's "no physics number defined
            # twice" applies to messages as well as constants.
            raise SystemExit(str(e))
        say("PYTHIA 8 bridge ready in %.1f s" % (time.time() - t0))
        _l.set_pythia_hadronizer(cfg, bridge)
        nthreads = 1                      # PythiaBridge is not re-entrant
    if opts["hfs_npz"] and not opts["hadronize"]:
        raise SystemExit("--hfs-npz needs --hadronize (no T2 hadrons "
                         "otherwise)")

    t0 = time.time()
    try:
        p = _l.Pipeline(cfg, plan)
    except RuntimeError as e:
        # CONFIGURATION-TIME REFUSALS ARE NOT TRACEBACKS.  Everything the
        # `Pipeline` constructor throws is a statement about the configuration
        # typed at this command line -- the run-plan spin against the ion
        # spin, a coherent window with no diffractive rate, and the one this
        # clause was written for: `InclusiveSampler`'s positivity check,
        # which evaluates 1 + w_avg over EVERY accepted cell before an event
        # is drawn and refuses when the phi-averaged density would go
        # negative.  It names the cell, the A1 = g1/F1 that made it negative
        # and the cure (`--x-max`); the `PythiaBridge` clause above is the
        # same rule, and the reason is the same one: the text is not
        # duplicated here (docs/CONVENTIONS.md, "no physics number is defined
        # twice", which applies to messages too).
        #
        # THE COMBINATION THAT FIRST NEEDED IT (measured 2026-09-05):
        # `--isotope d --channel tagged-d-p --unpol-sf ct18nlo --pol-sf
        # nnpdfpol --plan helicity-flip --pz 0.7 --pe 0.7` -- the A_par
        # measurement on the neutron target, at this CLI's own default --pe.
        # The neutron kernel's g1/F1 at the shipped window's top cell
        # (x = 0.955, Q2 = 1119) is 4.40579 on that backend PAIR against 0.7228
        # (ct18nlo + ToyG1) and 0.3820 (ToyF2 + nnpdfpol), so 1 + w_avg
        # = -0.02414 and the run is refused.  `--x-max 0.95` cures it.  It is
        # NOT clamped: a positivity clamp is a physics change.
        raise SystemExit(str(e))
    t_setup = time.time() - t0
    # THE RUN CONTEXT of the knob-provenance table: what the core cannot see.
    # The T2 fields are also recovered from the bound bridge inside the
    # metadata writer, which is the authority -- these are for the BANNER,
    # which is printed before a single event exists.
    ctx = _l.KnobRunContext()
    ctx.plan_name = opts["plan"]
    ctx.pz = float(opts["pz"])
    ctx.pzz = float(opts["pzz"])
    ctx.pzz_mode = str(opts["pzz_mode"])
    ctx.pe = float(opts["pe"])
    ctx.rel_lumi_offset = float(opts["rel_lumi_offset"])
    ctx.t2_bound = bridge is not None
    ctx.t2_pomeron = opts["coherent_t2"] != "off"
    ctx.pom_set = int(opts["pom_set"])
    ctx.pom_rescale = float(opts["pom_rescale"])
    # THE HEADER LINE GOES THROUGH THE TABLE TOO.  It printed
    # `p.optics.name` unconditionally until 2026-09-05, i.e. it announced the
    # far-forward envelope as a run property on EVERY channel -- including the
    # inclusive one, where the provenance block three lines below it says the
    # envelope is not read and the whole `route` column is Route::Lost.  One
    # banner cannot say both.  `KnobProvenance::meta_value` is the same rule
    # `meta["optics"]` now follows: the name where the envelope reached the
    # file, the scope clause where it did not.
    prov = p.knob_provenance(ctx)
    prov_by_name = {r.name: r for r in prov}
    say("%s  %s  %s" % (_l.pipeline_channel_name(cfg.channel),
                        p.beam_config.label(),
                        prov_by_name["optics"].meta_value))
    say("  setup %.2f s, %d events over %d categories, sigma = %.6g pb"
        % (t_setup, p.size(), len(plan), p.sigma_pb()))
    for name, sig, cnt in zip([c.name for c in plan.categories],
                              p.sigma_per_category_pb(), p.counts()):
        say("    %-10s sigma = %12.6g pb   N = %d" % (name, sig, cnt))
    # THE FILL'S OWN MOMENTS, on every run.  `--pzz` is read by the three
    # spin-1 tensor plans always and by `helicity-flip` only at
    # `--pzz-mode typed`: at the default `--pzz-mode ladder`
    # `helicity_flip_plan` leaves `use_explicit_pzz` false and takes the
    # MAX-ENTROPY ladder at `--pz`, so `--plan helicity-flip --pzz 0.6`
    # produces a different alignment from the one typed (T = 0.4 at pz = 0.7,
    # J = 3/2) and said so nowhere until 2026-09-05.  Printing what the plan
    # actually carries closes that without moving any fill: these are
    # `RunPlan`'s own recorded true moments, and at J = 3/2 `pzz_true` is the
    # normalised T, not P_zz (bookkeeping.hpp).
    #
    # THE SECOND LINE IS THE COUNTERFACTUAL, and it is DERIVED, never typed:
    # the ladder's own moment at this `--pz` comes from
    # `spin_temperature_pzz` (the same function `helicity_flip_plan`'s
    # default branch is built on) and the typed one is `--pzz` itself, so the
    # reader sees both fills and which one this run used.  The DOMAIN is not
    # printed here: it is stated once, in the refusal `spin1_populations` /
    # `spin32_populations` throw when it is crossed.
    j_fill = ion_spin(cfg.isotope)
    rank2 = "T" if abs(j_fill - 1.5) < 1e-9 else "P_zz"
    say("  fill %s: J = %g, P_z = %.6g, %s = %.6g"
        % (opts["plan"], j_fill, plan.pz_true, rank2, plan.pzz_true))
    if opts["plan"] in ("helicity-flip", "apar"):
        if opts["pzz_mode"] == "typed":
            say("    --pzz-mode typed: the fill is built at the TYPED "
                "%s = %.6g; the max-entropy ladder at this --pz would give "
                "%s = %.6g" % (rank2, opts["pzz"], rank2,
                               _l.spin_temperature_pzz(j_fill,
                                                       float(opts["pz"]))))
        else:
            say("    --pzz-mode ladder (the default): the max-entropy fill at "
                "--pz, so --pzz = %.6g is NOT read -- --pzz-mode typed would "
                "build the fill at %s = %.6g instead (refused, with the edge "
                "named, where that is outside the plan's domain)"
                % (opts["pzz"], rank2, opts["pzz"]))
    # THE COHERENT CHANNEL'S OWN LINES, on every coherent run.  The |t|
    # ceiling with the reason it is 0.2 (the anchor range, not positivity --
    # D5), the DERIVED positivity edge beside it so the two are visibly
    # different statements, and the model systematic of the T2 tier when the
    # tier is on.  All of it is in the npz `meta` too; the banner is where a
    # reader who never opens the file sees it.
    if cfg.channel == _l.PipelineChannel.CoherentLi6:
        sc = cfg.coherent
        # THE CONTINGENT NUMBER IS DERIVED, NOT TYPED.  Until 2026-09-05 this
        # line ended in the literal "2.80" -- no unit, no eps_b0 beside it,
        # and no basis -- while every other number on it was already derived.
        # It is also a rounding: 2.8000 is arithmetic on the ROUNDED
        # eps_b0 = -0.0070, and the DERIVED value is -0.0070024 -> 2.7990.
        # So build the measured-quadrupole scenario from the same band
        # `tests/test_coherent.cpp` T10b uses -- quadrupole_band_fm2()[0]
        # through a2_from_quadrupole -- at THIS run's own slope_b and amp,
        # and print what it returns.  (eps_b0 . B is the physical product, so
        # the derived eps_b0 scales with slope_b and the edge is invariant
        # under --coherent-slope-b; --coherent-amp does move it.)
        meas = _l.CoherentScenario()
        meas.slope_b = sc.slope_b
        meas.amp = sc.amp
        meas.eps_b0 = (_l.a2_from_quadrupole(
            2.0 * _l.ClusterConfigSampler().quadrupole_band_fm2()[0],
            6, 1.0, 1) * -4.0 / sc.slope_b)
        say("  coherent |t| <= %g GeV^2 (anchor range: the Mantysaari a2 "
            "table is digitized to |t| <= 0.30).  Positivity edge at "
            "eps_b0 = %g is |t| = %.4g (P_zz = -2) / %.4g (P_zz = +1) -- a "
            "CONTINGENT second reason: at the measured 6Li quadrupole "
            "(derived eps_b0 = %.5g) it is %.4f GeV^2"
            % (cfg.coherent_t_max, sc.eps_b0, sc.t_positivity_edge(-2.0),
               sc.t_positivity_edge(1.0), meas.eps_b0,
               meas.t_positivity_edge(-2.0)))
        if bridge is not None and opts["coherent_t2"] != "off":
            # WHETHER PomSet ran at all is the knob-provenance table's line
            # now (it was printed on this channel ONLY, so off it the banner
            # was silent while `meta` recorded a set that hadronized nothing).
            # What stays here is the BAND, which is a physics number and not a
            # reach claim.
            say("     T2 PomSet %d moves NOTHING at T0 (t, x_pom, q2, x, "
                "weight are bit-identical across the 14 sets that run: 1-10, "
                "12-15, set 11 refused); the band is on "
                "the HADRONIC FINAL STATE, over the 12 DPDF fits (3-10, "
                "12-15) about set 6: <n_ch> -2.6 %% / +9.9 %%, <n_had> "
                "-2.6 %% / +10.2 %%, <p_T> -5.9 %% / +10.9 %%, KAON FRACTION "
                "-44 %% / +55 %% (factor 2.8).  Run the sets and quote the "
                "envelope; never one set alone."
                % int(opts["pom_set"]))
    # WHICH STRUCTURE FUNCTIONS THIS RUN WAS MADE ON, on every run and every
    # channel, the all-default one included: the two backends move the rate
    # by 20 % and g1n by a sign, and a run that does not say which ones it
    # used is not reproducible from its own log.  Suppressed only when a
    # caller-supplied kernel wins, where the flags are not what produced the
    # numbers (`validate()` refuses a non-toy selector there anyway).
    if cfg.kernel is None:
        q2_min, below = _l.unpol_sf_grid_report(p)
        for line in sf_banner_lines(_l.unpol_sf_name(cfg.unpol_sf),
                                    _l.pol_sf_name(cfg.pol_sf),
                                    q2_min, below,
                                    _l.pol_sf_is_read(cfg, plan)):
            say(line)
    # THE KNOB-PROVENANCE BLOCK, on every run and every channel: ONE table
    # (`Pipeline.knob_provenance`) saying what this run read and what it did
    # not, replacing the per-knob reach sentences that were true of the knobs
    # somebody had looked at and silent about the rest.  The npz `meta`
    # carries the same table, from the same call.
    for line in knob_provenance_lines(prov):
        say(line)
    # THE LOUD ZERO, before the b1 block and unconditional on the run that
    # has it: a 7Li inclusive run's whole tensor and cos 2phi sector is
    # exactly 0.0, and until 2026-09-04 nothing on the run surface said so
    # (`inclusive_rank2_is_empty`, pipeline.hpp).  It is printed on EVERY such
    # run and not only on one with a tensor plan, for the same reason the SF
    # block above is: the user who has to be told is the one who did not know
    # to ask.
    if _l.inclusive_rank2_is_empty(cfg):
        for line in rank2_zero_lines(_l.rank2_input_report(cfg, plan),
                                     cfg.isotope):
            say(line)
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
                "-- the kernel's own f2_source is --unpol-sf %s, so at "
                "--unpol-sf toy the SPIN-BLIND cell cross section does not "
                "move with this flag; the tensor split does)"
                % (_l.b1_unpol_name(cfg.b1_unpol),
                   _l.unpol_sf_name(cfg.unpol_sf)))
            for line in b1_gate_lines(_l.b1_unpol_name(cfg.b1_unpol)):
                say("     " + line)
            # WHICH R the two halves of the tensor weight were built from.
            # Printed on every li6-convolution run, "unset" included, for the
            # `--b1-unpol` reason: at r1998 the weight moves by up to +26 %
            # at x = 0.10 and the unpolarised rate by -0.6847 %, and nothing
            # else in the banner says so.
            say("     R = sigma_L/sigma_T %s (--r-source; %s)"
                % (_l.r_source_name(cfg.r_source),
                   "unset = both hooks null, each resolving to r_sigma_lt "
                   "independently -- numerator and denominator agree today "
                   "only because two defaults coincide"
                   if cfg.r_source == _l.RSource.Unset else
                   "sigma-lt = the value both null hooks already resolve to; "
                   "bit-identical to unset, NOT READ (see KNOB PROVENANCE)"
                   if cfg.r_source == _l.RSource.SigmaLt else
                   "ONE object in BOTH Li6ConvolutionOptions::r_func (the "
                   "K numerator) and the kernel's Options::r_func (the "
                   "D_phi denominator, and F_L, D(y), ToyG1 -- so the "
                   "unpolarised rate moves with it)"))
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
            if r.options.tail_model == _l.RcTailModel.PolradFull:
                say("     the tail is POLRAD Eq. (18) + Appendix B + "
                    "Eq. (A.4) (--rc-tail-model polrad-full): ONE exact "
                    "tau_A quadrature carrying all THREE peaks, with the "
                    "s-/p-peaks' OWN Eq. (A.4) tensor content -- the term "
                    "t-peak+ll can only bound.  ITS ABSOLUTE NORMALISATION "
                    "IS STILL NOT CHECKED AGAINST MO-TSAI OR ANY EXTERNAL "
                    "EXACT TAIL: no such number is in this tree.  What is "
                    "checked is POLRAD-INTERNAL (the x_A -> 0 reduction to "
                    "Eq. (38), unpolarised AND tensor, and the Q_N = 0 "
                    "Rosenbluth limit) and the leading-log fallback.  The "
                    "QUASI-ELASTIC tail is a sum over SPIN-1/2 nucleons and "
                    "is tensor-blind at all three peaks here too.")
            elif r.options.tail_model == _l.RcTailModel.TPeak:
                say("     the tail is the t-PEAK ONLY (POLRAD Eqs. 37-39): "
                    "its absolute normalisation is NOT validated against an "
                    "exact tail.  Measured (test_rc.cpp T8(c), T8(d)) against "
                    "the leading-log s-/p-peaks, for 6Li at Q^2 >= 20 GeV^2 "
                    "and y <= 0.9: it carries the EVENT-WEIGHTED mean tail to "
                    "0.61 % at P_z = 0 (0.62 % at the CLI default P_z = 0.7), "
                    "but PER CELL only to 0.55 % over 0.15 <= y <= "
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
                    "tail: the tensor part of those peaks is BOUNDED, NOT "
                    "COMPUTED.  This is the UPPER edge of a band whose lower "
                    "edge is --rc-tail-model t-peak; quote both.")
                if cfg.rc_options.sp_tensor_scale == 0.0:
                    say("     the s/p TENSOR fraction is NOT priced on this "
                        "run: rc_tail carries the leading-log s-/p-peaks in "
                        "the UNPOLARISED numerator only, so this edge LOWERS "
                        "the tensor fraction of the tail by growing its "
                        "denominator alone.  --rc-sp-tensor-scale prices that "
                        "omission with a BOUND THAT HAS NO DERIVATION -- and "
                        "MEASURED (6Li config 1, 2026-09-06) that bound is "
                        "EMPTY at the three standard points: at x = 0.01 / "
                        "0.10 / 0.30, Q^2 = 5 the coherent s/p vertex sits at "
                        "Q'^2 = 4.37 / 4.94 / 4.98 GeV^2, where 6Li's form "
                        "factor is 45-51 decades down, so scale 1 is BIT-"
                        "IDENTICAL there.  It bites on 77 of 3051 accepted "
                        "cells (x < 0.008, y > 0.37), reaching 583 % of the "
                        "band half-width at x = 4.2e-4, y = 0.969.  The "
                        "QUASI-ELASTIC s/p column, which carries the s+p "
                        "everywhere else, is bounded by NOTHING.")
                else:
                    say("     the s/p TENSOR fraction is priced by a STAND-IN "
                        "at sp_tensor_scale = %.4g: the ELASTIC s-/p-peaks "
                        "are given the elastic t-peak's own tensor fraction. "
                        " IT IS A BOUND WITH NO DERIVATION -- the same "
                        "coherent 6Li vertex, but the s/p one sits at "
                        "Q'^2 ~ Q^2 where the coherent form factor is DEAD, "
                        "so the borrowed ratio may be the WRONG REFERENCE at "
                        "high Q^2, and there the term it scales is ~0, so it "
                        "prices NOTHING exactly where the s/p peaks dominate. "
                        " The QUASI-ELASTIC s/p column is covered by neither "
                        "this nor --rc-qe-tensor-scale and stays exactly "
                        "tensor-blind.  Quote it as a price tag on an "
                        "omission, and read the magnitude, never the sign."
                        % cfg.rc_options.sp_tensor_scale)
        if r.applies:
            say("     the weights are on Event.rc_weights and NOT on "
                "Event.weight: a systematic variation and a background, not "
                "a correction.")
            say("     Gakh-Shekhovtsova hep-ph/0403262 has ZERO INSPIRE "
                "citations: band it (--rc-delta-low-x 0.19 / 0.30, "
                "--rc-fq-scale 0 / 1 / 2, --rc-qe-suppression 0 / 0.5 / 1), "
                "never quote one row alone.")
            say("     and 0.30 is that paper's panel value carried UPWARD in "
                "x: at x = 0.00966, the x nearest this anchor's 0.01, the "
                "panel READS 0.113 (0.266 at its own bottom, x = 0.00226). "
                "The default errs WIDE by x2.65; priced at "
                "--rc-delta-low-x 0.266 / 0.113 in the 2026-09-06 run, "
                "sec. A2.")
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
    cols = p.generate(0, need_events, nthreads, context=ctx)
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
        # The Pomeron flavour fallback, beside n_ok/n_failed: it was counted
        # and surfaced NOWHERE (D4.7), which on a set that never consults the
        # Pomeron PDF is the difference between a run that used one and a run
        # that used none.  ~20 % is routine on the default set and costs
        # exactly nothing (pythia_bridge.hpp, PythiaBridgeStats).
        if st.n_pomeron:
            say("     Pomeron: %d hadronized, %d (%.2f %%) took the light "
                "e_q^2 flavour fallback -- measured bit-identical to the "
                "true draw (every Pomeron DPDF carries one light singlet), "
                "not a systematic"
                % (st.n_pomeron, st.n_pom_flavour_fallback,
                   100.0 * st.n_pom_flavour_fallback / st.n_pomeron))

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
        # Availability was checked before the run (see main()'s tier block).
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
