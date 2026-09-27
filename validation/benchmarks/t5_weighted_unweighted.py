#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""T5 MC closure -- ONE physics point generated UNWEIGHTED and in the tree's
WEIGHTED mode (Mode W), and the two cross-section estimates compared with
their MC errors (tier T5; the "partial" row of 05_chain.md sec. 1.1: "no
single test that generates the same physics point unweighted and weighted and
compares the two sigma estimates with their MC errors").

THE PHYSICS POINT.  The shipped defaults of `lipolgen.make_config()` --
6Li, inclusive, beam config 1 = `default_configs("6Li")[1]` (10 GeV e x
99.5 GeV/u), its generator scenario (Q2 >= 0.7, y <= 0.985, W2 >= 8,
|eta| <= 3.8, E' >= 0.3), a 100 x 72 log-log grid up to x = 1,
Q2 in [1, 2000], and `default_inclusive_kernel(li6())` (Miller b1, ToyF2,
ToyG1), which is what `Pipeline` builds for that config
(src/core/pipeline.cpp, inclusive branch) -- with TWO changes, both stated:
  * the window is the valence, high-depolarization corner x >= 0.1
    (`SamplerGrid.x_min`) and y >= 0.5 (`Scenario.y_min`);
  * the spin category is "apar+" of `make_plan("helicity-flip")` (J = 1,
    P_z = 0.7 ladder fill, lam_e = +1, P_e = 0.7).
WHY THE WINDOW.  On the full default window the polarized rate shift is
sigma(apar+)/sigma(unpol) - 1 = -6.80e-4, and a sigma comparison between an
unweighted and a weighted sample resolves it at 3 sigma only with >~ 4e7
events per route -- a closure that could not tell a working Mode W from
w = 1.  In this window the EXACT tables (sigma_tot_pb, noise-free) give a
shift of +1.47 %, resolved at 3 sigma by ~ 8.3e4 events per route; the
window was picked from those exact numbers alone (x >= 0.05/0.1/0.2 and
y >= 0.3/0.5/0.7 were tabulated), before any pull in it was computed.

THE TWO ROUTES, one fixed-seed stream each (SEED, runs 1 and 2, bunch 0):
  U  UNWEIGHTED (Mode G): `InclusiveSampler.sample_lumi(apar+, L)` -- the
     spin-labelled sampler draws the category's own density (m by
     p_m sigma_m, cells by sigma_c (1 + w_avg(m)), phi by accept-reject);
  W  WEIGHTED (Mode W, sampler.hpp "Mode-W style weights"): the SAME point
     generated UNPOLARIZED, `sample_lumi(unpol, L)` with populations
     (1/3, 1/3, 1/3), each event carrying w_i = `weights_for(ev, [apar+])`,
     the polarized/unpolarized density ratio at its cell and phi.

THE ESTIMATES AND THE TOLERANCE (the task's, fixed before any pull was
computed):
    sigma_U = N_U / L       +- sqrt(N_U) / L          (Poisson count)
    sigma_W = sum_i w_i / L +- sqrt(sum_i w_i^2) / L  (compound Poisson)
    pull = (sigma_U - sigma_W) / sqrt(err_U^2 + err_W^2),  PASS iff |pull| <= 3.
The two estimates share no event: U's count comes from the polarized
sampler's own normalisation (its m and cell CDFs), W's from the unpolarized
sampler times Mode W's weights.  They DO share the per-cell amplitude tables
(`StateTables` w_avg, a1, a2): an error IN those amplitudes moves both routes
alike and is not caught here -- this closes the sampling and the weighting
machinery, not the kernel's numbers.  Recorded, not asserted: each
estimate's pull against the exact `sigma_tot_pb(apar+)`; the fixed-N Mode-W
mean sigma_unpol <w> (what tests/test_sampler.cpp "Mode-W weights reproduce
the polarized rate ratio" checks) against the same; and the NULL pull a
no-op Mode W (w = 1) would give -- the closure's power.

N AND RUNTIME.  L = LUMI_PB = 34000 pb^-1 per nucleon per route (the window's
sigma_unpol is 87.9 pb), i.e. <N> ~ 2.99e6 events per route; single thread.
The pull's resolution is ~ 8.2e-4 relative, against the +1.47 % shift.
Measured 2.7 s wall for the whole harness, interpreter start-up included
(2.6 s inside run()), on the development container, 2026-09-26.

WHAT WAS MEASURED (2026-09-26, PASS).  Unweighted N_U/L = 89.1934 +- 0.0512
pb (N_U = 3 032 577) against weighted sum w/L = 89.1674 +- 0.0516 pb
(N_W = 2 987 706, <w> = 1.014722, w in [1.0061, 1.0676]): pull +0.359.
Exact sigma_tot_pb: 89.167123 pb (apar+), 87.873873 pb (unpol), shift
+1.4717e-2.  A no-op Mode W would pull +18.29.  Recorded: U vs exact +0.51,
W vs exact +0.00; the fixed-N Mode-W mean sigma_unpol <w> = 89.16753 +-
0.00043 pb, +0.94 against exact.

REFERENCE.  None vendored: an MC closure's reference is the tree's own second
route; no file under data/ is read.

THE THREE NORMALISATION RULES (BENCHMARK_PLAN.md sec. 8), for THIS row:
  1. b1 normalisations -- NOT APPLICABLE: no b1 is compared (the kernel's
     shipped Miller b1 is inside both routes' density identically).
  2. isoscalar denominator -- NOT APPLICABLE: no structure-function ratio.
  3. per nucleus / per nucleon -- both estimates are PER-NUCLEON cross
     sections in pb (xsec.hpp: F2A/A, per-nucleon x and s; the luminosity
     is per nucleon, L here in pb^-1 per nucleon), from ONE sampler, so the
     convention is common to both sides and cancels in the pull; nothing is
     converted.

WHAT THIS DOES NOT VALIDATE: the per-cell amplitudes themselves (shared
by both routes, above); the phi' part of Mode W (for this
longitudinal category a1 = a2 = 0, so phi is uniform and `weights_for`'s
cos phi', cos 2phi' terms are not exercised; the tensor a2 of the transverse
plans is ~ 1e-3); the rest of the phase space (the window is 1.5e-4 of the
default window's cross section; there the polarized shift is below any
affordable resolution); whether sigma_tot_pb is the true integral of the
cross section (the cell-centre grid discretisation is common to both routes);
the Pipeline's Event records (the harness drives the same InclusiveSampler
directly); the tagged, coherent and T1/T2 channels; the RC and FSI weight
columns; any physics of the polarized asymmetry.

HONESTY NOTE.  The tolerance |pull| <= 3 is the task's.  A FIRST VERSION of
this harness ran on the full default window (L = 5 pb^-1, ~2.96e6 events per
route) and measured count pull +0.365 (PASS), plus a second "shape" pair,
sigma_unpol/<1/w>_U against sigma_unpol <w>_W, pull -0.185 (PASS).  It was
replaced BEFORE this version's pull was computed, because (i) on that window
a no-op Mode W would have pulled only -0.46, and (ii) the shape pair was then
found, analytically, to be only SECOND-ORDER sensitive: its difference is
Cov(delta - eps, eps) ~ Var(w) ~ 1e-6 (delta, eps the sampler's and the
weights' polarized modulation), so it could not see even a sampler that
ignores the polarization.  The window above was chosen from the exact tables
only; the seed, the runs and the category did not change, and L was
rescaled to keep ~ 3e6 events per route.  The new pull (+0.359) is close to
the first version's (+0.365) NOT by coincidence of physics: both routes'
Poisson counts come from the same reserved counter of the same (seed, run,
bunch) streams, i.e. the same Poisson quantiles (U sits at +0.51 sigma of its
mean in both versions), so the two runs are not independent confirmations.
(The very first run died on a harness bug -- a duplicated `status` key --
after computing, before printing; the fix touched no number.)

Run:  source env.sh && python3 validation/benchmarks/t5_weighted_unweighted.py

Exit status (validation/benchmarks/README.md): 0 when the harness ran and
printed its REPORT row -- pass, recorded fail and blocked alike, since the
row carries the verdict; 2 when the harness itself is broken (a missing
dependency, any exception).  A harness is a measurement, not a CI gate: the
pytests are the gate.
"""

import math
import sys
import time
import traceback

NAME = "t5_weighted_unweighted"

SEED = 20260926
RUN_UNWEIGHTED = 1
RUN_WEIGHTED = 2
BUNCH = 0
X_MIN = 0.1                   # SamplerGrid.x_min of the window
Y_MIN = 0.5                   # Scenario.y_min of the window
LUMI_PB = 34000.0             # pb^-1 per nucleon, per route
PULL_MAX = 3.0
PLAN = "helicity-flip"
CATEGORY = "apar+"


def physics_point(x_min=X_MIN, y_min=Y_MIN):
    """The sampler on the window, the category, the unpolarized category."""
    import lipolgen as lg
    cfg = lg.make_config()
    beam = lg.default_configs(cfg.isotope)[cfg.beam_config]
    kern = lg.default_inclusive_kernel(lg.ion_by_name(cfg.isotope))
    scenario, grid = cfg.scenario, cfg.grid
    scenario.y_min = max(scenario.y_min, y_min)
    grid.x_min = max(grid.x_min, x_min)
    sampler = lg.InclusiveSampler(kern, beam, scenario, grid)
    cats = [c for c in lg.make_plan(PLAN).categories if c.name == CATEGORY]
    if len(cats) != 1:
        raise RuntimeError("%s: plan %r has no category %r" % (NAME, PLAN, CATEGORY))
    unpol = lg.SpinCategory("unpol", 1.0, [1.0 / 3.0] * 3)
    return dict(sampler=sampler, category=cats[0], unpol=unpol,
                isotope=cfg.isotope, beam=str(beam), beam_config=cfg.beam_config,
                x_min=float(grid.x_min), y_min=float(scenario.y_min))


def generate(pt, lumi_pb=LUMI_PB, seed=SEED):
    """Route U's count and route W's Mode-W weights."""
    s, cat, unpol = pt["sampler"], pt["category"], pt["unpol"]
    ev_u = s.sample_lumi(cat, lumi_pb, seed, run=RUN_UNWEIGHTED, bunch=BUNCH)
    n_u = int(ev_u["x"].size)
    del ev_u
    ev_w = s.sample_lumi(unpol, lumi_pb, seed, run=RUN_WEIGHTED, bunch=BUNCH)
    w_on_w = s.weights_for(ev_w, [cat])[:, 0]
    return n_u, w_on_w


def _pull(a, ea, b, eb):
    return (a - b) / math.sqrt(ea * ea + eb * eb)


def compare(n_u, w_on_w, lumi_pb, sigma_unpol=None, sigma_exact=None):
    """The two sigma estimates and their pull (pure numpy, so a test can hand
    it a mutated w)."""
    import numpy as np
    w = np.asarray(w_on_w, dtype=float)
    n_w = int(w.size)
    su, esu = n_u / lumi_pb, math.sqrt(n_u) / lumi_pb
    sw = float(w.sum()) / lumi_pb
    esw = math.sqrt(float((w * w).sum())) / lumi_pb
    out = dict(
        n_unweighted=int(n_u), n_weighted=n_w,
        sigma_u=su, err_u=esu, sigma_w=sw, err_w=esw,
        pull=_pull(su, esu, sw, esw),
        null_pull=_pull(su, esu, n_w / lumi_pb, math.sqrt(n_w) / lumi_pb),
        w_mean=float(w.mean()), w_std=float(w.std(ddof=1)),
        w_min=float(w.min()), w_max=float(w.max()),
    )
    if sigma_exact is not None:
        out["pull_u_vs_exact"] = (su - sigma_exact) / esu
        out["pull_w_vs_exact"] = (sw - sigma_exact) / esw
        if sigma_unpol is not None:
            mw = sigma_unpol * out["w_mean"]
            emw = sigma_unpol * out["w_std"] / math.sqrt(n_w)
            out["modew_mean"] = mw
            out["modew_mean_err"] = emw
            out["pull_modew_mean_vs_exact"] = (mw - sigma_exact) / emw
    out["status"] = "pass" if abs(out["pull"]) <= PULL_MAX else "fail"
    return out


def run(verbose=True, lumi_pb=LUMI_PB, seed=SEED, keep_weights=False):
    """The row.  `keep_weights` also returns route W's weight column (for a
    test that hands `compare` a mutated w); off, the report is numbers only."""
    t0 = time.perf_counter()
    pt = physics_point()
    s = pt["sampler"]
    sigma_unpol = float(s.sigma_tot_pb(pt["unpol"]))
    sigma_cat = float(s.sigma_tot_pb(pt["category"]))
    sigma_cells = float(sum(s.cell_xsec_pb))
    n_u, w_on_w = generate(pt, lumi_pb, seed)
    c = compare(n_u, w_on_w, lumi_pb, sigma_unpol, sigma_cat)
    runtime = time.perf_counter() - t0
    generator = ("unweighted sigma = %.4f +- %.4f pb (N = %d) vs weighted "
                 "(Mode W) %.4f +- %.4f pb (N = %d, <w> = %.5f): pull %+.3f"
                 % (c["sigma_u"], c["err_u"], c["n_unweighted"], c["sigma_w"],
                    c["err_w"], c["n_weighted"], c["w_mean"], c["pull"]))
    reference = ("same physics point, the other route: 6Li inclusive config %d, "
                 "x >= %g, y >= %g, helicity-flip %s (exact sigma_tot_pb %.4f pb; "
                 "a no-op Mode W would pull %+.1f)"
                 % (pt["beam_config"], pt["x_min"], pt["y_min"], CATEGORY,
                    sigma_cat, c["null_pull"]))
    tolerance = "|pull| <= %g" % PULL_MAX
    report = dict(
        name=NAME, generator=generator, reference=reference,
        tolerance=tolerance, status=c["status"],
        sigma_unpol_pb=sigma_unpol, sigma_cat_pb=sigma_cat,
        sigma_cells_pb=sigma_cells, rate_shift=sigma_cat / sigma_unpol - 1.0,
        lumi_pb=lumi_pb, seed=seed, runs=(RUN_UNWEIGHTED, RUN_WEIGHTED),
        category=CATEGORY, plan=PLAN, beam=pt["beam"], x_min=pt["x_min"],
        y_min=pt["y_min"], runtime_s=runtime,
        **{k: v for k, v in c.items() if k != "status"})
    if keep_weights:
        report["w_on_weighted"] = w_on_w
    if verbose:
        print("physics point: %s inclusive, %s, default kernel, generator scenario "
              "with y >= %g, grid x >= %g; category %s of %s"
              % (pt["isotope"], pt["beam"], pt["y_min"], pt["x_min"], CATEGORY, PLAN))
        print("exact: sigma_unpol = %.6f pb (sum of cells %.6f), sigma(%s) = %.6f pb, "
              "shift %+.4e" % (sigma_unpol, sigma_cells, CATEGORY, sigma_cat,
                               report["rate_shift"]))
        print("L = %g pb^-1/nucleon per route, seed %d, runs %d (U) / %d (W): "
              "N_U = %d, N_W = %d; w in [%.6f, %.6f], <w> = %.6f, std %.4e"
              % (lumi_pb, seed, RUN_UNWEIGHTED, RUN_WEIGHTED, c["n_unweighted"],
                 c["n_weighted"], c["w_min"], c["w_max"], c["w_mean"], c["w_std"]))
        print("  unweighted  N_U/L      = %.6f +- %.6f pb   (vs exact %+.2f)"
              % (c["sigma_u"], c["err_u"], c["pull_u_vs_exact"]))
        print("  weighted    sum w/L    = %.6f +- %.6f pb   (vs exact %+.2f)"
              % (c["sigma_w"], c["err_w"], c["pull_w_vs_exact"]))
        print("  pull (U - W) = %+.3f; a no-op Mode W (w = 1) would pull %+.2f"
              % (c["pull"], c["null_pull"]))
        print("  recorded: fixed-N Mode-W mean sigma_unpol <w> = %.6f +- %.6f pb "
              "(vs exact %+.2f)" % (c["modew_mean"], c["modew_mean_err"],
                                    c["pull_modew_mean_vs_exact"]))
        print("runtime %.2f s (in run())" % runtime)
        print("REPORT | %s | %s | %s | %s | %s" % (
            NAME, generator, reference, tolerance, c["status"].upper()))
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
