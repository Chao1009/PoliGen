# SPDX-License-Identifier: GPL-3.0-or-later
"""`lipolgen-configs` -- 6Li nucleon-position tables for a dipole-model code.

    lipolgen-configs --m +1 --n 100000 --theta-s 1.5707963 --phi-s 0 \
        --seed 20260902 --run 0 --out li6_m+1.dat \
        [--moments-json li6_m+1.moments.json]

Writes one line per configuration in the layout
`subnucleondiffraction`'s `Nucleons::InitializeTarget` reads for 3He
(positions in fm, ion rest frame, c.m. at the origin, one configuration per
line), with the polarization axis already applied, plus the mandatory sidecar
`<out>.meta.json`.  The C++ writer emits everything the core can compute; this
module fills in the two fields the core has no business implementing:

  * `md5` of every input table (hashlib), and
  * `git` -- `git rev-parse HEAD` when the enclosing work tree is LiPolGen's
    own checkout (verified through its `pyproject.toml`, so a venv sitting
    inside somebody else's repository does not get stamped in), else null.

With `--format annotated` the same JSON is ALSO copied into the `.dat` as a
`#` header block by the C++ writer -- that copy is written before this module
runs, so its `md5`/`git` stay null.  The sidecar is the authoritative record.

It is also the ONLY place a human-facing warning can live: `src/core` has no
`std::cerr` and setup errors throw (docs/CONVENTIONS.md), so `--theta-s 0`
prints "a_2 == 0 for a longitudinal axis" from here while C++ simply returns
0, which is the correct value there.

CAVEAT, PRINTED ON EVERY RUN.  The alpha+d truncation reproduces the 6Li point
radius to ~4 % but OVERSHOOTS Q(6Li): with the DEFAULT source Q_charge(model)
= -0.615 fm^2 against the measured -0.0818 and GFMC AV18+IL7's -0.20(6), a
factor ~7.5.  The printed factor is COMPUTED from the run's own
`quadrupole_band_fm2()`, so it follows `--alpha-d-source`, `--alpha-d-scale`
and `--quadrupole-target` instead of restating one number; the three sources'
own values are design_G_cluster_config.md sec. 8.  Never quote a tensor number
from these tables without the band; see
docs/open_items/run_2026-09-02/design_G_cluster_config.md sec. 2.7.

That ~7.5 is TWO factors, not three (measured 2026-09-03,
docs/open_items/run_2026-09-03/phase_C_numbers.md sec. C1): 3.3165 from the
model to the dial setting that matches the MEASURED asymptotic D/S ratio, times
2.2686 from there to the measurement.  1/S_alpha-d = 1.1706 is NOT a third
factor -- both waves are divided by sqrt(S_alpha-d) before the moments are
taken, so it is already inside the -0.615.  And the leg's error bar is the
physics: George & Knutson's +-0.011662 on eta maps, through a dial that eta is
EXACTLY linear in, to a model Q from -0.4005 fm^2 to +0.0298 fm^2 -- through
ZERO.  Quote it as 3.32x (1 sigma: 1.54x .. sign change) or not at all.

AND `a2_at_t0p3` IS A COEFFICIENT, NOT A SENSITIVITY (measured 2026-09-04,
phase_C_numbers.md sec. C2 / OPEN_ITEMS_SOLUTIONS.md sec. 11.3; arithmetic in
validation/o5_a2_reach.py).  It is a_2 at |t| = 0.3 GeV^2 and the coherent
sample lives at |t| ~ 1/B ~ 0.026, so the modulation the EIC would see is
kappa*sqrt(<t^2>) = 0.25 % at the MEASURED Q, a factor 10.6 smaller.  Priced:
coherent J/psi over the whole Q^2 range with both lepton channels is MARGINAL
at one EIC year, inside the {1, 10, 100} fb^-1/u band.  QUOTE IT AS A BAND:
S = 2.63 sigma at the band's LOW EDGE and 2.84 .. 3.29 sigma at its TOP, with
3 sigma at 8.3 .. 13.0 fb^-1/u (sec. 11.3b).  The detection efficiency in that
chain is a 7Li number at 7Li's own TOP energy used at 10 x 99.5 (x1.12-1.16
UP, on the same paper's own 3He energy scan) while the chain carries NO
decay-lepton acceptance or reconstruction efficiency at all (DOWN, unbounded
below in this tree).  Those two cancel to 0.7 %, which is why the uncorrected
point estimate -- 2.62 sigma, 3 sigma at 13.1 fb^-1/u -- lands 0.3 % under the
band's low edge, and why neither the point nor either edge may be quoted
alone.  WHETHER THE BAND'S TOP CROSSES 3 SIGMA IS NOT ESTABLISHED: the top is
a span because the 7Li -> 6Li efficiency substitution straddles 1 once it is
read off entries that share a beam energy (sec. 11.3b, `o5_a2_reach.py`).
(An earlier revision of this docstring said "NOT MEASURABLE ... 0.75 sigma,
160 fb^-1/u"; that was one lepton channel in one Q^2 window,
OPEN_ITEMS_SOLUTIONS.md sec. 11.3a.)

READ THAT NUMBER WITH ITS LIMITATION.  It comes from `a2_from_quadrupole`,
which is a CLOSED FORM AND NOT A GOOD-WALKER DIPOLE-MODEL AMPLITUDE: the
target's quadrupole through the deuteron's published |t| dependence, with no
amplitude, no saturation and none of their uncertainties, and with the MATTER
quadrupole standing in for the transverse GLUON anisotropy.  FOUR things are
also unestablished rather than uncertain: no detection efficiency exists below
Q^2 = 0.1 GeV^2 anywhere in this tree, where 85 % of the rate sits; no
decay-lepton reconstruction efficiency exists in this tree at all (only its
geometric half is bounded, at 0.99), which is what leaves the band open below;
the 7Li -> 6Li efficiency substitution straddles 1, which is what makes the
band's top a span; and the far-forward working point is unchosen -- the LAST
is the single correction that on its own restores the NO (at LiPolGen's own
de-squeezed 6Li tagging optics the band is 0.73-0.92 sigma with 3 sigma at
106-167 fb^-1/u).

`lipolgen-run` is untouched: this is a separate console entry point, not a
flag on it, because adding subparsers there would change its parse behaviour.
"""

import argparse
import hashlib
import json
import math
import os
import subprocess
import sys
import time

from . import _lipolgen as _l

ALPHA_SOURCES = {"vmc": _l.AlphaCoreSource.VmcHe4Density,
                 "gaussian": _l.AlphaCoreSource.Gaussian}
ALPHA_D_SOURCES = {"fit-rescaled": _l.AlphaDSource.FitRescaled,
                   "fit-raw": _l.AlphaDSource.FitRaw,
                   "overlap-raw": _l.AlphaDSource.OverlapRaw}
FORMATS = {"he3": _l.SndConfigFormat.He3Compatible,
           "annotated": _l.SndConfigFormat.Annotated}

DEFAULTS = dict(m="+1", n=100000, theta_s=math.pi / 2.0, phi_s=0.0,
                seed=20260902, run=0, alpha_source="vmc",
                alpha_d_source="fit-rescaled", alpha_d_scale=1.0,
                match_li6_radius=False, quadrupole_target=0.0,
                min_nn_separation=0.0, alpha_cm_inflate=True,
                n_r=512, n_c=96, r_max=20.0, rnp_max=25.0,
                fmt="he3", out=None, moments_json=None, quiet=False)


def build_parser():
    p = argparse.ArgumentParser(
        prog="lipolgen-configs",
        description="Write 6Li nucleon-position configuration tables for a "
                    "Good-Walker dipole-model code (subnucleondiffraction).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    p.add_argument("--m", default=DEFAULTS["m"],
                   choices=["+1", "1", "0", "-1", "unpolarized"],
                   help="ion substate; 'unpolarized' interleaves equal thirds")
    p.add_argument("--n", type=int, default=DEFAULTS["n"],
                   help="number of configurations")
    p.add_argument("--theta-s", type=float, default=DEFAULTS["theta_s"],
                   help="quantization-axis polar angle [rad]; a cos 2Phi "
                        "signal needs theta_s != 0")
    p.add_argument("--phi-s", type=float, default=DEFAULTS["phi_s"])
    p.add_argument("--seed", type=int, default=DEFAULTS["seed"])
    p.add_argument("--run", type=int, default=DEFAULTS["run"])
    p.add_argument("--alpha-source", default=DEFAULTS["alpha_source"],
                   choices=sorted(ALPHA_SOURCES))
    p.add_argument("--alpha-d-source", default=DEFAULTS["alpha_d_source"],
                   choices=sorted(ALPHA_D_SOURCES))
    p.add_argument("--alpha-d-scale", type=float,
                   default=DEFAULTS["alpha_d_scale"])
    p.add_argument("--match-li6-radius", action="store_true",
                   help="override --alpha-d-scale with the value that puts "
                        "<r^2> on li6.density's 2.4433 fm")
    p.add_argument("--quadrupole-target", type=float,
                   default=DEFAULTS["quadrupole_target"],
                   help="Q_charge [fm^2] to dial the alpha-d D wave to "
                        "((G9) root, NOT sqrt(target/model)); 0 = off")
    p.add_argument("--min-nn-separation", type=float,
                   default=DEFAULTS["min_nn_separation"],
                   help="hard core between the alpha's nucleons [fm]")
    p.add_argument("--no-alpha-cm-inflate", dest="alpha_cm_inflate",
                   action="store_false",
                   help="do not inflate the alpha source before recentring "
                        "(then <s^2> comes out 3/4 of the table's)")
    p.add_argument("--n-r", type=int, default=DEFAULTS["n_r"])
    p.add_argument("--n-c", type=int, default=DEFAULTS["n_c"])
    p.add_argument("--r-max", type=float, default=DEFAULTS["r_max"])
    p.add_argument("--rnp-max", type=float, default=DEFAULTS["rnp_max"])
    p.add_argument("--format", dest="fmt", default=DEFAULTS["fmt"],
                   choices=sorted(FORMATS),
                   help="'he3' is byte-compatible with he3.dat (no comments); "
                        "'annotated' needs the patched upstream reader")
    p.add_argument("--out", default=None, help="configuration table to write")
    p.add_argument("--moments-json", default=None,
                   help="write the analytic and sampled moments side by side")
    p.add_argument("--quiet", action="store_true")
    p.set_defaults(alpha_cm_inflate=DEFAULTS["alpha_cm_inflate"])
    return p


def make_options(opts):
    """A `ClusterConfigOptions` from the resolved option dict."""
    o = _l.ClusterConfigOptions()
    o.alpha_source = ALPHA_SOURCES[opts["alpha_source"]]
    o.alpha_d_source = ALPHA_D_SOURCES[opts["alpha_d_source"]]
    o.theta_s = float(opts["theta_s"])
    o.phi_s = float(opts["phi_s"])
    o.alpha_cm_inflate = bool(opts["alpha_cm_inflate"])
    o.min_nn_separation_fm = float(opts["min_nn_separation"])
    o.alpha_d_scale = float(opts["alpha_d_scale"])
    o.quadrupole_target_fm2 = float(opts["quadrupole_target"])
    o.n_r = int(opts["n_r"])
    o.n_c = int(opts["n_c"])
    o.r_max_fm = float(opts["r_max"])
    o.rnp_max_fm = float(opts["rnp_max"])
    return o


def _git(here, *args):
    """`git -C here <args>` -> stripped stdout, or None."""
    try:
        out = subprocess.run(["git", "-C", here, *args],
                             capture_output=True, text=True, timeout=10)
    except (OSError, subprocess.SubprocessError):
        return None
    return out.stdout.strip() or None if out.returncode == 0 else None


def _is_lipolgen_tree(top):
    """Does `top` look like LiPolGen's own checkout?"""
    try:
        with open(os.path.join(top, "pyproject.toml")) as f:
            head = f.read(8192)
    except OSError:
        return False
    return 'name = "lipolgen"' in head


def git_head(path=None):
    """`git rev-parse HEAD` of LiPolGen's OWN checkout, else None.

    The enclosing repository is CHECKED, not assumed.  A non-editable wheel
    puts this file in site-packages, and a `.venv` commonly lives inside the
    user's own analysis repository -- a bare `git -C <here> rev-parse HEAD`
    would then stamp THAT repository's SHA into the sidecar as the
    generator's provenance.  A SHA is returned only when the enclosing
    work tree's `pyproject.toml` is LiPolGen's own (which is still true for
    the staged `build/python/lipolgen` copy, whose toplevel is this repo).
    """
    here = path or os.path.dirname(os.path.abspath(__file__))
    top = _git(here, "rev-parse", "--show-toplevel")
    if top is None or not _is_lipolgen_tree(top):
        return None
    return _git(here, "rev-parse", "HEAD")


def file_md5(path):
    h = hashlib.md5()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def fill_sidecar(meta_path):
    """Fill the `md5` and `git` fields the C++ writer emits as null."""
    with open(meta_path) as f:
        meta = json.load(f)
    meta["git"] = git_head()
    for entry in meta.get("inputs", []):
        path = _l.data_path(entry["file"])
        try:
            entry["md5"] = file_md5(path)
        except OSError:
            entry["md5"] = None
    with open(meta_path, "w") as f:
        json.dump(meta, f, indent=1)
        f.write("\n")
    return meta


def moments(sampler, cfg_set, elapsed_s=None):
    """The analytic predictors and the sampled moments, side by side."""
    n = len(cfg_set)
    band = sampler.quadrupole_band_fm2()
    out = dict(
        analytic=dict(
            p_d_alpha_d=sampler.p_d_alpha_d(),
            s_alpha_d=sampler.s_alpha_d(),
            tensor_dilution=sampler.tensor_dilution(),
            r2_alpha_fm2=sampler.r2_alpha_fm2(),
            r2_np_fm2=sampler.r2_np_fm2(),
            r2_alpha_d_fm2=sampler.r2_alpha_d_fm2(),
            q_int_fm2=sampler.q_int_fm2(),
            q_dd_fm2=sampler.q_dd_fm2(),
            quadrupole_fm2=sampler.q_int_fm2() + sampler.q_dd_fm2(),
            r2_mean_fm2=sampler.r2_analytic_fm2(),
            r_rms_fm=math.sqrt(sampler.r2_analytic_fm2()),
            q_matter_fm2={str(m): sampler.q_matter_analytic_fm2(m)
                          for m in (1, 0, -1)},
            q_charge_fm2=0.5 * sampler.q_matter_analytic_fm2(1),
            delta_perp_fm2=sampler.delta_perp_analytic_fm2(1),
            eps_b0_equivalent=sampler.eps_b0_equivalent(),
            a2_at_t0p3={str(m): sampler.a2_from_geometry(0.3, m)
                        for m in (1, 0)},
            asymptotic_ds_ratio=sampler.asymptotic_ds_ratio(),
            quadrupole_dial_s=sampler.quadrupole_dial_s(),
            match_li6_radius_scale=sampler.match_li6_radius(),
        ),
        grid=dict(r2_mean_fm2=sampler.r2_grid_fm2(),
                  q_matter_fm2=sampler.q_matter_grid_fm2(1),
                  p2_alpha_d_m1_ms_minus1=sampler.p2_alpha_d_grid(1, -1.0),
                  p2_alpha_d_m1_ms_zero=sampler.p2_alpha_d_grid(1, 0.0)),
        sampled=dict(n=n, m_ion=cfg_set.m_ion, seed=cfg_set.seed,
                     run=cfg_set.run,
                     r2_mean_fm2=cfg_set.r2_mean_fm2,
                     r_rms_fm=math.sqrt(cfg_set.r2_mean_fm2) if n else 0.0,
                     r2_sd_fm2=cfg_set.r2_sd_fm2,
                     q_matter_fm2=cfg_set.q_matter_fm2,
                     q_matter_sd_fm2=cfg_set.q_matter_sd_fm2,
                     delta_perp_fm2=cfg_set.delta_perp_fm2,
                     delta_perp_sd_fm2=cfg_set.delta_perp_sd_fm2,
                     cm_max_fm=cfg_set.cm_max),
        quadrupole_band_fm2=dict(measured=band[0], gfmc_av18_il7=band[1],
                                 gfmc_av18_il7_err=_l.LI6_QUADRUPOLE_GFMC_ERR_FM2,
                                 alpha_d_model=band[2]),
    )
    if n > 1:
        s = math.sqrt(float(n))
        out["sampled"]["five_sigma"] = dict(
            r2_mean_fm2=5.0 * cfg_set.r2_sd_fm2 / s,
            q_matter_fm2=5.0 * cfg_set.q_matter_sd_fm2 / s,
            delta_perp_fm2=5.0 * cfg_set.delta_perp_sd_fm2 / s)
    if elapsed_s is not None:
        out["timing"] = dict(wall_s=elapsed_s,
                             us_per_config=1e6 * elapsed_s / n if n else 0.0)
    return out


def resolve(argv=None):
    args = vars(build_parser().parse_args(argv))
    opts = dict(DEFAULTS)
    opts.update({k: v for k, v in args.items() if v is not None})
    return opts


def main(argv=None):
    opts = resolve(argv)
    say = (lambda *a: None) if opts["quiet"] else (lambda *a: print(*a))

    if float(opts["theta_s"]) == 0.0:
        say("lipolgen-configs: WARNING -- a_2 == 0 for a longitudinal axis "
            "(theta_s = 0): the projected density is azimuthally symmetric, "
            "so there is no cos 2Phi modulation to measure.  Use "
            "--theta-s 1.5707963 for a transverse axis.")

    o = make_options(opts)
    if opts["match_li6_radius"]:
        probe = _l.ClusterConfigSampler(o)
        o.alpha_d_scale = probe.match_li6_radius()
        say("lipolgen-configs: --match-li6-radius -> alpha_d_scale = %.6f"
            % o.alpha_d_scale)
    sampler = _l.ClusterConfigSampler(o)

    band = sampler.quadrupole_band_fm2()
    say("lipolgen-configs: CAVEAT -- Q_charge(model) = %+.4f fm^2 vs measured "
        "%+.4f and GFMC AV18+IL7 %+.2f +- %.2f fm^2."
        % (band[2], band[0], band[1], _l.LI6_QUADRUPOLE_GFMC_ERR_FM2))
    say("                  |model/measured| = %.1f; carry the band, never one "
        "number." % abs(band[2] / band[0]))
    if sampler.quadrupole_dial_s() != 1.0 or o.alpha_d_scale != 1.0:
        say("                  DIALLED/SCALED (dial s = %.4f, alpha_d_scale = "
            "%.4f): that factor is this TUNED geometry's, not the natural "
            "wave functions' overshoot."
            % (sampler.quadrupole_dial_s(), o.alpha_d_scale))
    if o.min_nn_separation_fm > 0.0:
        say("lipolgen-configs: NOTE -- min_nn_separation %.2f fm: the alpha's "
            "<s^2> is a NUMERICAL estimate (the closed form (3/4)<v^2> holds "
            "only for independent draws), so the analytic <r^2> / r_rms and "
            "match_li6_radius() are approximate and the sidecar marks them."
            % o.min_nn_separation_fm)

    t0 = time.time()
    if opts["m"] == "unpolarized":
        cfg_set = sampler.sample_set_unpolarized(int(opts["n"]),
                                                 int(opts["seed"]),
                                                 int(opts["run"]))
    else:
        m = int(opts["m"].lstrip("+"))
        cfg_set = sampler.sample_set(int(opts["n"]), m, int(opts["seed"]),
                                     int(opts["run"]))
    elapsed = time.time() - t0
    say("lipolgen-configs: %d configurations in %.3f s (%.2f us each)"
        % (len(cfg_set), elapsed,
           1e6 * elapsed / max(len(cfg_set), 1)))

    if opts["out"]:
        rows = _l.write_snd_configs(cfg_set, sampler, opts["out"],
                                    FORMATS[opts["fmt"]])
        meta_path = opts["out"] + ".meta.json"
        fill_sidecar(meta_path)
        say("lipolgen-configs: wrote %d rows to %s (+ %s)"
            % (rows, opts["out"], meta_path))

    mom = moments(sampler, cfg_set, elapsed)
    if opts["moments_json"]:
        with open(opts["moments_json"], "w") as f:
            json.dump(mom, f, indent=1)
            f.write("\n")
        say("lipolgen-configs: moments -> %s" % opts["moments_json"])
    elif not opts["out"]:
        say(json.dumps(mom, indent=1))
    return 0


if __name__ == "__main__":
    sys.exit(main())
