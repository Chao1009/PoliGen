# SPDX-License-Identifier: GPL-3.0-or-later
"""LiPolGen -- doubly polarized e + 6Li / 7Li DIS event generator.

The C++ core (`_lipolgen`, pybind11) is re-exported wholesale, so everything
in `docs/USAGE.md` is spelled the same way here:

    import lipolgen as lg

    cfg = lg.PipelineConfig()
    cfg.channel  = lg.PipelineChannel.TaggedLi6Alpha
    cfg.isotope  = "6Li"
    cfg.beam_config = 1                      # 10 GeV e x 99.5 GeV/u
    cfg.n_events = 100000
    cfg.seed     = 1
    p = lg.Pipeline(cfg, lg.tensor_thirds_plan(0.7, 0.6))

    ev = p.generate()                        # columnar dict of numpy arrays
    print(ev["x"].mean(), (ev["route"] == 1).mean())

`Pipeline.generate(n=0, events=False, nthreads=1)` returns a dict of numpy
arrays -- one row per event, every array owning its own buffer -- plus a
`meta` dict.  `Pipeline.generate_lumi(lumi_pb)` rebuilds the same run at an
integrated luminosity and generates all of it.

The convenience layer on top (`make_config`, `make_plan`, `run`) is what the
`lipolgen-run` command line uses, and is the shortest path from a few strings
to a sample:

    cols = lg.run(isotope="6Li", channel="tagged-alpha", plan="tensor-thirds",
                  events=100000, seed=1)

`lipolgen.export` turns either form into the polligen schemas
(`inclusive_dict`, `tagged_dict`, `write_hfs_npz`).
"""

# If a wheel-packaged copy of data/vmc lives next to this file (see the
# `install(DIRECTORY ... DESTINATION lipolgen/data)` rule in CMakeLists.txt)
# and the caller has not already pointed $LIPOLGEN_DATA_DIR somewhere,
# prefer it.  `_lipolgen.data_dir()` reads the environment once, lazily, on
# first use (src/core/cluster.cpp), so this only has to run before that --
# doing it here, before the extension is even imported, is early enough.
# In-tree builds (`build/python/lipolgen/`) never have a `data/` sibling, so
# this is a no-op there and the compiled-in `$CMAKE_SOURCE_DIR/data` default
# keeps resolving exactly as before.
import os as _os

if not _os.environ.get("LIPOLGEN_DATA_DIR"):
    _pkg_data = _os.path.join(_os.path.dirname(__file__), "data")
    if _os.path.isdir(_pkg_data):
        _os.environ["LIPOLGEN_DATA_DIR"] = _pkg_data
    del _pkg_data

from ._lipolgen import *          # noqa: F401,F403
from . import _lipolgen           # noqa: F401
from ._lipolgen import __version__

from . import export              # noqa: F401
from .export import (inclusive_dict, tagged_dict, hfs_sample,  # noqa: F401
                     write_hfs_npz, write_columns_npz)

__all__ = [n for n in dir(_lipolgen) if not n.startswith("_")] + [
    "export", "inclusive_dict", "tagged_dict", "hfs_sample", "write_hfs_npz",
    "write_columns_npz", "CHANNELS", "PLANS", "TRITON_SFS", "FSI", "RC",
    "B1_MODELS", "B1_UNPOL", "R_SOURCE", "UNPOL_SF", "POL_SF", "PZZ_MODES",
    "RC_C0_SHAPES", "RC_TAIL_MODELS", "OPTICS", "CLUSTER_WAVES",
    "make_config", "make_plan", "make_pipeline", "run", "ion_spin",
    "__version__",
]
# `OPTICS`, `CLUSTER_WAVES` and `ion_spin` were missing from the list until
# 2026-09-16 while every sibling name table was on it and `cli.py` imports all
# three, so `from lipolgen import *` gave FSI and RC but not OPTICS.


# --------------------------------------------------------------- name tables

#: Command-line channel names -> (PipelineChannel factory, allowed isotopes).
#: "tagged-alpha" resolves against the isotope, exactly as
#: `pipeline.channel_isotope` does in reverse.
CHANNELS = {
    "inclusive": lambda iso: _lipolgen.PipelineChannel.Inclusive,
    "tagged-alpha": lambda iso: (_lipolgen.PipelineChannel.TaggedLi6Alpha
                                 if iso == "6Li"
                                 else _lipolgen.PipelineChannel.TaggedLi7Alpha),
    "tagged-6Li-alpha": lambda iso: _lipolgen.PipelineChannel.TaggedLi6Alpha,
    "tagged-7Li-alpha": lambda iso: _lipolgen.PipelineChannel.TaggedLi7Alpha,
    "tagged-d-p": lambda iso: _lipolgen.PipelineChannel.TaggedDeuteronP,
    "coherent": lambda iso: _lipolgen.PipelineChannel.CoherentLi6,
}

#: Command-line optics names -> OpticsChoice.
OPTICS = {
    "yr-high-acceptance": _lipolgen.OpticsChoice.YellowReportHighAcceptance,
    "yr-high-divergence": _lipolgen.OpticsChoice.YellowReportHighDivergence,
    "tagging": _lipolgen.OpticsChoice.Tagging,
    "tagging-legacy": _lipolgen.OpticsChoice.TaggingLegacyLevers,
}

#: cluster radial-form families.  "hulthen" is the DEFAULT and is
#: bit-compatible with every published number; "vmc" swaps in the ANL VMC
#: tables (see docs/CONVENTIONS.md and docs/open_items/vmc_reconciliation.md)
#: and then ignores `cluster_beta` and `p_d` on the lithium alpha tags.
#: It is read on the TAGGED channels only and `PipelineConfig.validate`
#: refuses it elsewhere; on the 6Li alpha tag it selects the embedded
#: deuteron as well as the alpha-d relative motion (C5.5b).
CLUSTER_WAVES = {
    "hulthen": _lipolgen.ClusterWaveSource.Hulthen,
    "vmc": _lipolgen.ClusterWaveSource.VmcAV18,
}

#: triton spectral functions of the 7Li alpha tag's T1 breakup.  "hulthen" is
#: the DEFAULT (the sequential two-body decay, bit-compatible with every
#: published 7Li number); "ciofi-simula" swaps in the Ciofi degli Atti-Simula
#: n0 + n1 spectral function (`triton_sf.hpp`, docs/CONVENTIONS.md): the
#: k-dependent n0/(n0+n1) branching integrating to S0 = 0.6525 and the third
#: channel (struck n -> a (p n) continuum) the sequential model does not have.
TRITON_SFS = {
    "hulthen": _lipolgen.TritonSfChoice.Hulthen,
    "ciofi-simula": _lipolgen.TritonSfChoice.CiofiSimula,
}

#: FSI weight models accepted on the command line.  "off" (the default) is
#: today's plane-wave impulse approximation bit for bit; the other two apply
#: the Glauber FSI/IA WEIGHT to Event.weight (never a shift of any
#: four-vector; spin independent by construction -- see fsi.hpp).  Quote the
#: effect as an unpolarized-shape systematic, never as a correction to A_zz,
#: and BAND sigma_XN over 20-40 mb rather than quoting one row alone.
FSI = {
    "off": _lipolgen.PipelineFsi.Off,
    "glauber-cluster": _lipolgen.PipelineFsi.GlauberCluster,
    "glauber-nucleon": _lipolgen.PipelineFsi.GlauberNucleon,
}

#: RC weight families accepted on the command line.  "off" (the default) is
#: today bit for bit; "tensor-band" emits rc_tensor_lo / rc_tensor_hi /
#: rc_tail on Event.rc_weights and on NOTHING else -- never Event.weight,
#: never a four-vector, never a random number.  The band is a SYSTEMATIC
#: VARIATION and the tail a BACKGROUND, so an analysis multiplies them in on
#: purpose; band the knobs (--rc-delta-low-x 0.19 / 0.30, --rc-fq-scale
#: 0 / 1 / 2, --rc-qe-suppression 0 / 0.5 / 1) and never quote one row alone.
RC = {
    "off": _lipolgen.RcMode.Off,
    "tensor-band": _lipolgen.RcMode.TensorBand,
}

#: b1 backends of the INCLUSIVE kernel's rank-2 slot (`--b1-model`).
#: "miller" is the DEFAULT and is bit-for-bit what every published inclusive
#: tensor number was made with -- Li6B1(MillerB1) through
#: LI6_B1_RANK2_TRANSFER -- and it IS the run PLAN's "toy": today's default
#: b1_func is Li6B1(MillerB1), which is what `toy_b1` reaches.  "cdks" is the
#: same 6Li rank-2 transfer on the other CAMP for b1_d (the digitized CDKS
#: PRD 95 (2017) 074036 Fig. 4 column: |b1| two orders of magnitude smaller
#: below x ~ 0.1, comparable or larger above it -- peak |x b1| 3.34e-4 at
#: x = 0.766 against Miller's 4.27e-4 at x = 0.084, both at Q2 = 2.5, and 7.7x
#: LARGER at x = 0.3 with the opposite sign) --
#: Miller and CDKS disagree and the library does not adjudicate, so say which
#: one a plot used.  INCLUSIVE CHANNEL ONLY and 6Li ONLY, like
#: "li6-convolution".  "li6-convolution" is
#: b1_nuclear.hpp's four-term alpha-d convolution (design_D_b1_li6.md):
#: INCLUSIVE CHANNEL ONLY and 6Li ONLY (on a tagged channel the alpha-d
#: density is already in the event weight).  Its A = 2 magnitude gate PASSES
#: since 2026-09-03, but only for ONE unpolarised nucleon input: MSTW2008 LO
#: (`--b1-unpol mstw`, `B1_UNPOL` below) at CDKS Eq. (21)'s delta-function,
#: where G3b = 0.843243.  The SHIPPED DEFAULT `--b1-unpol toy` gives 0.440 --
#: OUTSIDE the gate's [0.5, 2] acceptance window -- so quote numbers made with
#: `b1_unpol="mstw"`; toy-backend numbers are not covered by the lift.  The
#: gate is A = 2 on every setting and says nothing about the alpha-d step, so
#: band every number with --b1-band-scale 0/1/2 and never quote one row alone.
B1_MODELS = {
    "miller": _lipolgen.B1Model.Miller,
    "cdks": _lipolgen.B1Model.Cdks,
    "li6-convolution": _lipolgen.B1Model.Li6Convolution,
}

#: `--b1-unpol` -- which UNPOLARISED structure-function backend the
#: "li6-convolution" b1 folds its own F1 against
#: (`Li6ConvolutionOptions.unpol`).  "toy" is the DEFAULT and is bit for bit
#: what every published number was made with: the inclusive kernel's own
#: `ToyF2`, handed to the convolution as ONE shared object.  "mstw" is
#: MSTW2008 LO over PYTHIA 8's own pdfdata grid -- the PDF CDKS computed their
#: b1_d with, and the one the A = 2 gate passes on (G3b 0.843 against 0.440 on
#: the toy) -- and needs the optional PYTHIA tier.  "ct18nlo" is the phase-D
#: stand-in, kept selectable so the systematic can be quoted rather than
#: remembered, and needs the optional LHAPDF tier.
#:
#: MEASURED on the shipped observable `Li6ConvolutionB1.b1(x, 2.5)`:
#: mstw/toy = 1.847766 / 1.275961 / 0.816971 and ct18nlo/toy = 2.221302 /
#: 1.238833 / 1.045870 at x = 0.10 / 0.30 / 0.50.  Up to a factor 1.85 and NOT
#: monotone, so it is a shape change and not a normalisation.
#:
#: It reaches the b1 and NOTHING else.  The kernel's own `f2_source` is set
#: by a SEPARATE flag, `--unpol-sf` (`UNPOL_SF` below), and at its shipped
#: default "toy" it is `ToyF2` on every `--b1-unpol` setting -- which is what
#: keeps the SPIN-BLIND cell cross section (`InclusiveSampler.cell_xsec_pb`)
#: bit for bit under THIS flag, with only the tensor shift moving.  With
#: `unpol_sf` set to anything else the rate moves too, but it moves because
#: of that flag; `validate()` refuses `b1_unpol="toy"` there, because "toy"
#: means "the kernel's own UnpolSF, shared" and that object would no longer
#: be `ToyF2`.  Read only by "li6-convolution"; `validate()` refuses it
#: on the other two models rather than let `meta["b1_unpol"]` record a PDF
#: that never touched the rate.  `B1UnpolSource.Custom` is not here on
#: purpose: it names an object, so it is set by assigning
#: `config.b1_unpol_sf`.
B1_UNPOL = {
    "toy": _lipolgen.B1UnpolSource.Toy,
    "mstw": _lipolgen.B1UnpolSource.Mstw,
    "ct18nlo": _lipolgen.B1UnpolSource.Ct18Nlo,
}

#: `--r-source` -- the ONE R = sigma_L/sigma_T threaded into BOTH halves of
#: the 6Li tensor weight.  The shipped observable is a RATIO:
#: A_zz = -(2/3) K/D_phi with K = b1 + (1-y)/(x y^2) b2 and
#: D_phi = F1 + (1-y)/(x y^2) F2 (`asymmetries.azz`).  Its numerator's R is
#: `Li6ConvolutionOptions.r_func` and its denominator's is the inclusive
#: kernel's own `Options.r_func`; both default to null, which `resolve_r`
#: turns into `r_sigma_lt`, so they agree today by coincidence of two
#: independent defaults and the only way to move one was to move it ALONE.
#: This selector is registry option (iii) of `STATUS.md` row 3: whichever R
#: is named, `default_inclusive_kernel` puts the SAME object in both.
#:
#: "unset" is the DEFAULT and installs NOTHING -- both hooks stay null and
#: the run is bit for bit the pre-flag tree, by construction.  "sigma-lt"
#: installs one `r_sigma_lt` in both and is MEASURED bit-identical to
#: "unset" (sigma_pb, all three per-category cross sections and every
#: generated column; and b1, K/D_phi, A_zz and the cos 2phi amplitude at the
#: six standard points x = 0.05/0.10/0.30 x Q2 = 2.5/5, each at y = 0.1, 0.5
#: and 0.9): it is the wiring's own test, not a physics variation, and
#: `knob_provenance` labels it `not-read` for exactly that reason.  "r1998" installs the A = 2 gate's R in both.
#:
#: MEASURED at y = 0.5, Q2 = 2.5 (docs/open_items/run_2026-09-06/
#: phase_A_numbers.md sec. A1), against "unset": the tensor weight K/D_phi --
#: and A_zz, which is -(2/3) of it -- moves -3.84 % / +26.40 % / +3.35 % at
#: x = 0.05 / 0.10 / 0.30, and the cos 2phi amplitude, whose numerator
#: carries no R at all, moves -8.33 % / -6.73 % / -3.14 % there.  Unlike
#: `--b1-unpol` this one moves the UNPOLARISED rate too (sigma_pb -0.6847 %
#: on a 2000-event 6Li run at x_max 0.95), because the kernel's `r_func`
#: reaches F1, F_L, D(y) and `ToyG1`.
#:
#: Read only by "li6-convolution"; `validate()` refuses it on the other two
#: models rather than let `meta["r_source"]` record an R that never touched
#: the rate.
R_SOURCE = {
    "unset": _lipolgen.RSource.Unset,
    "sigma-lt": _lipolgen.RSource.SigmaLt,
    "r1998": _lipolgen.RSource.R1998,
}

#: `--unpol-sf` -- which UNPOLARISED structure-function backend supplies F2,
#: and through it F1, F_L and the whole unpolarised rate, to EVERY kernel a
#: `Pipeline` builds: the inclusive kernel, the coherent channel that rides
#: its cell cross sections, and the tagged struck-cluster kernel.  It is the
#: injection point the tagged channels did not have at all before 2026-09-04
#: -- `PipelineConfig.kernel` is read on the non-tagged branch alone.
#:
#: "toy" is the DEFAULT and is bit for bit what every published number was
#: made with: `ToyF2`, which `sf.hpp` labels TOY and anchors BY EYE
#: ("adequate for phase-space maps and factor-1.5 rate estimates ONLY").
#: "mstw" is MSTW2008 LO over PYTHIA 8's own pdfdata grid (optional PYTHIA
#: tier); "ct18nlo" is `LhapdfSF("CT18NLO", 0)` (optional LHAPDF tier).
#:
#: MEASURED, 6Li inclusive at config 1 over the shipped window: the summed
#: accepted cell cross section is 591846.2 pb (toy), 472571.9 (ct18nlo,
#: x0.7985), 469556.5 (mstw, x0.7934).  It is a shape change, not a
#: normalisation: F2p at Q2 = 10 moves x0.9215 / x1.1095 / x1.3709 / x1.2641
#: (ct18nlo) and x0.8154 / x1.0133 / x1.3649 / x1.4479 (mstw) at
#: x = 0.01 / 0.10 / 0.30 / 0.50, and the toy's F2n/F2p -- the straight line
#: clip(1 - 0.75x, 0.25, 1), which the T1 and T2 species draws use -- is
#: 0.9625 / 0.8500 / 0.6250 at x = 0.05 / 0.20 / 0.50 against CT18NLO's
#: 0.9218 / 0.7219 / 0.5035.
#:
#: READ THE GRID CLAUSE for "ct18nlo": CT18NLO's grid starts at Q2 = 1.677
#: and the shipped window's accepted cells start at 1.054, so 36.18 % of a
#: CT18NLO run's own accepted cell cross section (42.33 % of the toy run's)
#: sits below the grid, where LHAPDF keeps evolving downward instead of
#: freezing.  It is not refused -- that would make the flag unusable on the
#: shipped scenario -- but it is printed at the banner and recorded in
#: `meta["unpol_sf_below_grid_frac"]` on every run.
#:
#: NOT ORTHOGONAL TO `POL_SF`: `InclusiveKernel` builds its default `ToyG1`
#: on its own base `UnpolSF`, so this flag alone moves g1 as well and
#: A1 = g1/F1 moves 5-8 %.  `pol_sf="toy"` does not mean "g1 unchanged".
#: `UnpolSfSource.Custom` is not here on purpose: it names an object, so it
#: is set by assigning `config.unpol_sf_obj`.
UNPOL_SF = {
    "toy": _lipolgen.UnpolSfSource.Toy,
    "mstw": _lipolgen.UnpolSfSource.Mstw,
    "ct18nlo": _lipolgen.UnpolSfSource.Ct18Nlo,
}

#: `--pol-sf` -- which POLARISED backend supplies g1, and through the
#: Wandzura-Wilczek relation g2, to the INCLUSIVE and TAGGED kernels a
#: `Pipeline` builds.  NOT to every kernel, and unlike `UNPOL_SF` it never
#: could: g1 enters through the single product lam_e * P_e * (m/J) *
#: cos(theta_S) * A_par, so the selector is read only where the CHANNEL
#: evaluates g1 AND the FILL carries lam_e * P_e != 0.  The coherent channel
#: does neither, and neither does any unpolarised-beam plan -- "tensor-thirds"
#: is the CLI's own default and builds every category at lam_e = 0, so at the
#: defaults `--pol-sf` is read on NO channel.  It is LABELLED there, not
#: refused and not credited (`pol_sf_is_read(config, plan)`, and the `pol_sf`
#: row of `Pipeline.knob_provenance`).  "toy"
#: is the DEFAULT and is bit for bit: `ToyG1` on the kernel's own `UnpolSF`
#: and `r_func`.  "nnpdfpol" is `LhapdfG1("NNPDFpol11_100", 0)` -- already
#: `LhapdfG1`'s own declared default, and the ONLY polarised set installed in
#: this tree's LHAPDF store, so there is no second row to offer.
#:
#: MEASURED on 6Li per-nucleon g1A (A_par tracks it to better than 0.1 % at
#: y = 0.5): nnpdfpol/toy = 0.6691 / 0.6761 / 1.1139 / 1.3491 / 1.2982 /
#: 0.9815 / 0.6358 at (x, Q2) = (0.01, 2.5) / (0.05, 5) / (0.10, 10) /
#: (0.20, 10) / (0.30, 15) / (0.50, 25) / (0.70, 50) -- x0.64 to x1.35 over
#: the generator window, and not monotone.
#:
#: THE NEUTRON IS A SIGN, NOT A FACTOR.  `ToyG1`'s
#: a1n(x) = -0.07(1-x)^2 + 0.8 x^2.2 crosses zero near x ~ 0.25 and is
#: POSITIVE above it, while NNPDFpol1.1's g1n stays negative to x ~ 0.6:
#: g1n(x, 10) is -0.0801 / -0.0114 / +0.00521 / +0.00779 (toy) against
#: -0.1365 / -0.0608 / -0.02740 / -0.00037 (nnpdfpol) at
#: x = 0.10 / 0.20 / 0.30 / 0.50.  The shipped toy g1n has the WRONG SIGN
#: over roughly 0.25 < x < 0.6.  On isoscalar 6Li the proton term dominates
#: and this mostly hides; on a neutron-tagged run (d + p tagging, whose
#: `dis_target` is NEUTRON_TARGET) it does not.
POL_SF = {
    "toy": _lipolgen.PolSfSource.Toy,
    "nnpdfpol": _lipolgen.PolSfSource.NnpdfPol,
}

#: The two edges of the 6Li C0 (monopole) shape band, for `rc_c0_shape` /
#: `--rc-c0-shape`.  It is shared by `F_c` and `F_q`, so it moves BOTH and
#: cannot be rescaled out of one run the way a multiplier can.
#:
#: "ho" is the DEFAULT and is bit for bit what every published number was made
#: with: the unfitted (a, alpha) harmonic oscillator, whose F_c and F_q both
#: change sign at q_0 = 3.0999 fm^-1.  "vmc-ft" is the j0 transform of the
#: committed ANL VMC point-proton density, r-rescaled to the SAME measured
#: <r^2>_point = 6.0788 fm^2 (so F_c(0), F_q(0) and <r^2> are identical on the
#: two edges) and exponentially continued above q = 3 fm^-1, so it has NO C0
#: zero at any q.
#:
#: MEASURED on `(1/6) sigma^el_T / sigma^el_U` at Q^2 = 5 GeV^2:
#: ho/vmc-ft = -5.0943e-4 / -5.4705e-4 at x = 0.01, +1.5595e-4 / -4.0702e-5 at
#: x = 0.10 and +1.3287e-2 / +1.5556e-2 at x = 0.30 -- i.e. the SIGN at
#: x = 0.10 is a band edge, not a result (phase_B_numbers.md sec. B1).
#: RUN BOTH.
RC_C0_SHAPES = {
    "ho": _lipolgen.C0Shape.Ho,
    "vmc-ft": _lipolgen.C0Shape.VmcFt,
}

#: The THREE tail formulations, for `rc_tail_model` / `--rc-tail-model`.
#: (`RcTailModel.PolradFull` was refused by `PipelineConfig.validate` and
#: absent from this table until 2026-09-06; it is implemented now.)
#:
#: "t-peak" is the DEFAULT and is bit for bit every published number: POLRAD
#: Eqs. (37)-(39), (43), the t-peak ALONE, and therefore a LOWER BOUND on the
#: dilution.  "t-peak+ll" adds the leading-log s- and p-peaks of the same two
#: unpolarised observables.
#:
#: IT IS A STATED MODEL, NOT A CONTROLLED EXPANSION.  The t-peak is an eta_A
#: quadrature; the s-/p-peaks are a single-z collinear leading log.  Their sum
#: double-counts nothing but is accurate only to the leading log, ~5-10 %.
#:
#: AND IT HAS NO TENSOR PARTNER.  POLRAD's Eq. (38) supplies no sigma_T at the
#: s-/p-peak and this library will not invent one from the leading log, so the
#: s+p contribution enters the UNPOLARISED numerator only.  (Eq. (18) +
#: Eq. (A.4) DOES supply one and "polrad-full" computes it: the unqualified
#: sentence is true of Eq. (38) and false of the paper, which is why
#: rc_sp_tensor_scale is refused on "polrad-full" because the term RAN.)
#: Switching to "t-peak+ll" therefore LOWERS the
#: tensor FRACTION of the tail (x0.66 / x0.0032 / x6.6e-05 at x = 0.01 / 0.10 /
#: 0.30, Q^2 = 5): BOUNDED by rc_sp_tensor_scale, not computed.  RUN BOTH.
#:
#: "polrad-full" is POLRAD Eq. (18) + Appendix B + Eq. (A.4): ONE exact tau_A
#: quadrature that contains all three peaks, with the s-/p-peaks' OWN tensor
#: content -- the thing "t-peak+ll" has to bound and cannot compute.  It is
#: NOT the default and NOT validated against any external exact tail (no
#: Mo-Tsai number is in this tree): what is checked is its x_A -> 0 reduction
#: to Eq. (38), unpolarised AND tensor, its Q_N = 0 Rosenbluth limit, and its
#: agreement with the leading-log fallback where the s-/p-peaks are alive.
#: It costs ~10 s of table build and needs n_eta >= 64 and a `long double`
#: wider than `double`; both are refused, not degraded.
RC_TAIL_MODELS = {
    "t-peak": _lipolgen.RcTailModel.TPeak,
    "t-peak+ll": _lipolgen.RcTailModel.TPeakPlusLL,
    "polrad-full": _lipolgen.RcTailModel.PolradFull,
}

#: Run-plan names accepted on the command line (aliases included).
PLANS = ("tensor-thirds", "azz", "helicity-flip", "apar", "transverse-tensor",
         "cos2phi", "tensor-flip", "flip")

#: `--pzz-mode` -- WHICH FILL the `helicity-flip` / `apar` plan builds, and
#: the only knob in this module that names a branch of a plan FACTORY rather
#: than a field of `PipelineConfig`.
#:
#:   "ladder"  the DEFAULT.  `helicity_flip_plan` leaves
#:             `HelicityFlipOptions.use_explicit_pzz` false and takes the
#:             SPIN-TEMPERATURE (max-entropy) populations at `pz`, so `pzz`
#:             is not read at all and the fill carries the rank-2 moment a
#:             pure vector fill of J >= 1 necessarily drags along
#:             (`spin_temperature_pzz`).  Bit for bit the tree before this
#:             switch existed, by construction: it is the branch that was
#:             already there.
#:   "typed"   `use_explicit_pzz = true`: the fill is
#:             `spin1_populations(pz, pzz)` at J = 1 and
#:             `spin32_populations(pz, pzz)` at J = 3/2, i.e. the typed
#:             alignment, with the octupole R_3 left at 0.  A value outside
#:             the plan's domain is REFUSED with the edge named, never
#:             clamped: 3|P_z| - 2 <= P_zz <= 1 at J = 1 and the smaller
#:             1.8|P_z| - 1 <= T <= 1 - 0.6|P_z| at J = 3/2, so at P_z = 0.7
#:             the typed T = 0.6 of the CLI default is outside the J = 3/2
#:             domain by 0.02 (the edge is T = 0.58, where p(-1/2) = 0).
#:
#: The two fills are DIFFERENT PHYSICS, not a correction of one another, and
#: the cost of choosing the second is measured at the standard configuration
#: in `docs/open_items/run_2026-09-06/phase_A_numbers.md` sec. A3.  The three
#: tensor plans have no ladder branch, so the mode leaves their fills exactly
#: as they are and `knob_provenance` labels it `not-read` there.
PZZ_MODES = ("ladder", "typed")

#: Ion spin J by species -- the spin the run plan has to be built for.
_ION_SPIN = {"p": 0.5, "d": 1.0, "3He": 0.5, "6Li": 1.0, "7Li": 1.5}


def ion_spin(isotope):
    """Ion spin J of a species name ("6Li" -> 1.0, "7Li" -> 1.5)."""
    return _ION_SPIN[isotope]


# ------------------------------------------------------------- construction

def make_config(isotope="6Li", config=1, channel="inclusive", events=0,
                lumi_pb=0.0, seed=20260713, run=1, optics="yr-high-acceptance",
                poisson=True, cluster_beta=None, p_d=None, coherent=None,
                scenario=None, grid=None, n_sigma=10.0, pot_config="",
                inclusive_b1=None, with_virtual_photon=True,
                apply_optics_lumi_fraction=None,
                cluster_wave=None, triton_sf=None,
                fsi=None, fsi_sigma_mb=None,
                rc=None, rc_delta_low_x=None, rc_delta_high_x=None,
                rc_a_transfer_frac=None,
                rc_fq_scale=None, rc_tail_tensor_scale=None,
                rc_qe_suppression=None, rc_qe_tensor_scale=None,
                rc_sp_tensor_scale=None,
                rc_c0_shape=None,
                rc_tail_model=None,
                b1_model=None, b1_band_scale=None,
                b1_alpha_d_dwave_weight=None, b1_unpol=None,
                r_source=None,
                unpol_sf=None, pol_sf=None, coherent_t_max=None):
    """A `PipelineConfig` from plain values (the CLI's own constructor).

    `channel` is a key of `CHANNELS`; `optics` a key of `OPTICS`.  The isotope
    a channel implies wins over the `isotope` argument (a 6Li alpha tag is a
    6Li run whatever the caller said), which is `channel_isotope`'s rule.
    `coherent` is a dict of `CoherentScenario` fields, and `coherent_t_max`
    the channel's |t| ceiling [GeV^2] (`PipelineConfig.coherent_t_max`, NOT a
    `CoherentScenario` field): 0.2 by default, kept there by the Mantysaari
    anchor's |t| <= 0.30 digitization range rather than by positivity, and
    recorded in the npz meta.  `cluster_wave` is a
    key of `CLUSTER_WAVES` ("hulthen", the default, or "vmc") or a
    `ClusterWaveSource` directly; "vmc" replaces the lithium alpha-tag radial
    forms with the ANL VMC tables and then ignores `cluster_beta` / `p_d`.
    `triton_sf` is a key of `TRITON_SFS` ("hulthen", the default, or
    "ciofi-simula") or a `TritonSfChoice` directly; "ciofi-simula" replaces
    the 7Li alpha tag's sequential triton breakup with the three-channel
    Ciofi degli Atti-Simula spectral function, built by the Pipeline at the
    run's own `cluster_beta`.
    `fsi` is a key of `FSI` ("off", the default, "glauber-cluster" or
    "glauber-nucleon") or a `PipelineFsi` directly, and `fsi_sigma_mb` the
    sigma_XN it is built at (40 = free hadron; band 20-40 mb -- run both
    ends, never quote one row alone).  Tagged channels only; the FSI enters
    as a per-event weight on `Event.weight` and moves no four-vector.
    `rc` is a key of `RC` ("off", the default, or "tensor-band") or an
    `RcMode` directly, and the `rc_*` knobs write the matching
    `RcOptions` fields -- including `rc_a_transfer_frac`, the A = 2 -> A = 6
    transfer price on the band (design Q8; 0.0 default = the shipped band).
    The RC weights land on `Event.rc_weights` and on
    NOTHING else -- never `Event.weight`, never a four-vector, never a random
    number -- so "off" is today bit for bit.
    `b1_model` is a key of `B1_MODELS` ("miller", the default, "cdks" or
    "li6-convolution") or a `B1Model` directly; it chooses what
    `default_inclusive_kernel` puts in the INCLUSIVE kernel's rank-2 slot and
    is ignored when a caller-supplied `kernel` is set (`validate()` refuses
    that combination for anything but "miller").  Both opt-in backends are
    INCLUSIVE CHANNEL ONLY and 6Li ONLY -- "cdks" is `Li6B1`'s 6Li rank-2
    transfer too -- and `validate()` refuses the rest rather than let the
    metadata record a flag that never reached the rate.  `b1_band_scale` is
    the MANDATORY 100 % band on those two backends -- run 0 / 1 / 2 and quote
    the envelope, never a single row -- and `b1_alpha_d_dwave_weight` the
    shape knob on the alpha-d orbital terms (2d) + (2a) together, read only by
    "li6-convolution".  A knob that the chosen backend does not read is
    REFUSED at 1.0-away values, for the same provenance reason.
    `r_source` is a key of `R_SOURCE` ("unset", the default, "sigma-lt" or
    "r1998"): the ONE R = sigma_L/sigma_T threaded into BOTH the alpha-d
    convolution's own F1 and the kernel's, so the tensor weight's numerator
    and denominator are the same choice.  Read by
    `b1_model="li6-convolution"` only.

    `b1_unpol` is a key of `B1_UNPOL` ("toy", the default, "mstw" or
    "ct18nlo") or a `B1UnpolSource` directly: the UNPOLARISED backend the
    "li6-convolution" b1 folds its own F1 against.  "mstw" is CDKS's own
    MSTW2008 LO and is the configuration the A = 2 gate passes on; it needs
    the optional PYTHIA tier and raises RuntimeError naming it (or naming the
    missing grid file) in a build that has none -- it is never silently
    replaced by the toy.  It is read only by "li6-convolution" and moves b1
    only: AT `unpol_sf="toy"` the spin-blind cell cross section is
    bit-identical across its settings, while the tensor-weighted per-category
    cross sections are not.
    `unpol_sf` is a key of `UNPOL_SF` ("toy", the default, "mstw" or
    "ct18nlo") and `pol_sf` a key of `POL_SF` ("toy", the default, or
    "nnpdfpol"), or the matching enum directly: the structure-function
    backends of the kernels the run builds.  THEY DO NOT HAVE THE SAME REACH.
    `unpol_sf` reaches EVERY kernel -- inclusive, coherent and tagged -- on
    every plan.  `pol_sf` reaches the inclusive and tagged kernels only where
    the fill also carries lam_e * P_e != 0, and the coherent channel not at
    all: under "tensor-thirds", this function's own default plan, it is read
    on no channel and the run records the LABEL rather than the backend name
    (`pol_sf_is_read`).  Both default to the toy backends and are then bit for
    bit.  A
    named backend needs its optional tier and raises RuntimeError naming the
    missing one; it is never silently replaced by the toy, which is x0.80 on
    the 6Li rate and has the WRONG SIGN on g1n over 0.25 < x < 0.6.  They are
    NOT orthogonal -- `unpol_sf` alone moves g1 through the default `ToyG1`,
    so state both next to any number.  `unpol_sf` is a separate choice from
    `b1_unpol` (that one is the deuteron F1 inside CDKS Eq. (22), a
    CDKS-comparability choice); `validate()` refuses the one combination that
    would mislabel them, `b1_unpol="toy"` with a non-toy `unpol_sf` on
    "li6-convolution".
    """
    if channel not in CHANNELS:
        raise ValueError("unknown channel %r; know %s"
                         % (channel, ", ".join(sorted(CHANNELS))))
    cfg = _lipolgen.PipelineConfig()
    cfg.channel = CHANNELS[channel](isotope)
    implied = _lipolgen.channel_isotope(cfg.channel)
    cfg.isotope = implied or isotope
    cfg.beam_config = int(config)
    cfg.seed = int(seed)
    cfg.run = int(run)
    cfg.poisson = bool(poisson)
    cfg.n_sigma = float(n_sigma)
    cfg.pot_config = pot_config
    cfg.with_virtual_photon = bool(with_virtual_photon)
    if events:
        cfg.n_events = int(events)
    if lumi_pb:
        cfg.lumi_pb = float(lumi_pb)
    if optics is not None:
        if optics not in OPTICS:
            raise ValueError("unknown optics %r; know %s"
                             % (optics, ", ".join(sorted(OPTICS))))
        cfg.optics_choice = OPTICS[optics]
    if cluster_beta is not None:
        cfg.cluster_beta = float(cluster_beta)
    if p_d is not None:
        cfg.p_d = float(p_d)
    if cluster_wave is not None:
        cfg.cluster_wave = CLUSTER_WAVES[cluster_wave] \
            if cluster_wave in CLUSTER_WAVES else cluster_wave
    if triton_sf is not None:
        if isinstance(triton_sf, str):
            if triton_sf not in TRITON_SFS:
                raise ValueError("unknown triton_sf %r; know %s"
                                 % (triton_sf, ", ".join(sorted(TRITON_SFS))))
            cfg.triton_sf = TRITON_SFS[triton_sf]
        else:
            cfg.triton_sf = triton_sf
    if fsi is not None:
        if isinstance(fsi, str):
            if fsi not in FSI:
                raise ValueError("unknown fsi %r; know %s"
                                 % (fsi, ", ".join(sorted(FSI))))
            cfg.fsi = FSI[fsi]
        else:
            cfg.fsi = fsi
    if fsi_sigma_mb is not None:
        cfg.fsi_sigma_mb = float(fsi_sigma_mb)
    if rc is not None:
        if isinstance(rc, str):
            if rc not in RC:
                raise ValueError("unknown rc %r; know %s"
                                 % (rc, ", ".join(sorted(RC))))
            cfg.rc = RC[rc]
        else:
            cfg.rc = rc
    # `RcOptions` is returned BY VALUE from the binding, so every knob has to
    # be written on one copy and assigned back -- the `cfg.struck` pattern.
    if any(v is not None for v in (rc_delta_low_x, rc_delta_high_x,
                                   rc_a_transfer_frac,
                                   rc_fq_scale, rc_tail_tensor_scale,
                                   rc_qe_suppression, rc_qe_tensor_scale,
                                   rc_sp_tensor_scale,
                                   rc_c0_shape, rc_tail_model)):
        opt = cfg.rc_options
        if rc_delta_low_x is not None:
            opt.delta_low_x = float(rc_delta_low_x)
        if rc_delta_high_x is not None:
            opt.delta_high_x = float(rc_delta_high_x)
        if rc_a_transfer_frac is not None:
            opt.a_transfer_frac = float(rc_a_transfer_frac)
        if rc_fq_scale is not None:
            opt.fq_scale = float(rc_fq_scale)
        if rc_tail_tensor_scale is not None:
            opt.tail_tensor_scale = float(rc_tail_tensor_scale)
        if rc_qe_suppression is not None:
            opt.qe_suppression = float(rc_qe_suppression)
        if rc_qe_tensor_scale is not None:
            opt.qe_tensor_scale = float(rc_qe_tensor_scale)
        if rc_sp_tensor_scale is not None:
            opt.sp_tensor_scale = float(rc_sp_tensor_scale)
        if rc_c0_shape is not None:
            if isinstance(rc_c0_shape, str):
                if rc_c0_shape not in RC_C0_SHAPES:
                    raise ValueError("unknown rc_c0_shape %r; know %s"
                                     % (rc_c0_shape,
                                        ", ".join(sorted(RC_C0_SHAPES))))
                opt.c0_shape = RC_C0_SHAPES[rc_c0_shape]
            else:
                opt.c0_shape = rc_c0_shape
        if rc_tail_model is not None:
            if isinstance(rc_tail_model, str):
                if rc_tail_model not in RC_TAIL_MODELS:
                    raise ValueError("unknown rc_tail_model %r; know %s"
                                     % (rc_tail_model,
                                        ", ".join(sorted(RC_TAIL_MODELS))))
                opt.tail_model = RC_TAIL_MODELS[rc_tail_model]
            else:
                opt.tail_model = rc_tail_model
        cfg.rc_options = opt
    if scenario is not None:
        cfg.scenario = scenario
    if grid is not None:
        cfg.grid = grid
    if apply_optics_lumi_fraction is not None:
        cfg.apply_optics_lumi_fraction = bool(apply_optics_lumi_fraction)
    if b1_model is not None:
        if isinstance(b1_model, str):
            if b1_model not in B1_MODELS:
                raise ValueError("unknown b1_model %r; know %s"
                                 % (b1_model, ", ".join(sorted(B1_MODELS))))
            cfg.b1_model = B1_MODELS[b1_model]
        else:
            cfg.b1_model = b1_model
    if b1_band_scale is not None:
        cfg.b1_band_scale = float(b1_band_scale)
    if b1_alpha_d_dwave_weight is not None:
        cfg.b1_alpha_d_dwave_weight = float(b1_alpha_d_dwave_weight)
    if b1_unpol is not None:
        if isinstance(b1_unpol, str):
            if b1_unpol not in B1_UNPOL:
                raise ValueError("unknown b1_unpol %r; know %s"
                                 % (b1_unpol, ", ".join(sorted(B1_UNPOL))))
            b1_unpol = B1_UNPOL[b1_unpol]
        # NOT a plain assignment: the enum is the provenance and the backend
        # object is the realisation, and the core library links neither the
        # PYTHIA nor the LHAPDF tier, so the binding is what builds the
        # object.  It raises here, naming the missing tier, rather than
        # letting `validate()` report a named-but-empty slot later.
        _lipolgen.set_b1_unpol(cfg, b1_unpol)
    if r_source is not None:
        if isinstance(r_source, str):
            if r_source not in R_SOURCE:
                raise ValueError("unknown r_source %r; know %s"
                                 % (r_source, ", ".join(sorted(R_SOURCE))))
            r_source = R_SOURCE[r_source]
        # A PLAIN ASSIGNMENT, unlike the three selectors around it: both R
        # functions live in the CORE (`sf.hpp`), so there is no optional tier
        # to build and no provenance/realisation pair to keep together.
        cfg.r_source = r_source
    if unpol_sf is not None:
        if isinstance(unpol_sf, str):
            if unpol_sf not in UNPOL_SF:
                raise ValueError("unknown unpol_sf %r; know %s"
                                 % (unpol_sf, ", ".join(sorted(UNPOL_SF))))
            unpol_sf = UNPOL_SF[unpol_sf]
        # NOT a plain assignment, for the `set_b1_unpol` reason above: the
        # enum is the provenance and the object the realisation, and only
        # this layer can see the optional tiers.
        _lipolgen.set_unpol_sf(cfg, unpol_sf)
    if pol_sf is not None:
        if isinstance(pol_sf, str):
            if pol_sf not in POL_SF:
                raise ValueError("unknown pol_sf %r; know %s"
                                 % (pol_sf, ", ".join(sorted(POL_SF))))
            pol_sf = POL_SF[pol_sf]
        _lipolgen.set_pol_sf(cfg, pol_sf)
    if inclusive_b1 is not None:
        struck = cfg.struck
        struck.inclusive_b1 = bool(inclusive_b1)
        cfg.struck = struck
    if coherent:
        sc = cfg.coherent
        for k, v in coherent.items():
            if not hasattr(sc, k):
                raise ValueError("CoherentScenario has no field %r" % k)
            setattr(sc, k, float(v))
        cfg.coherent = sc
    if coherent_t_max is not None:
        # NOT a `CoherentScenario` field -- it is `PipelineConfig`'s own |t|
        # ceiling (coherent.hpp, COHERENT_T_MAX_DEFAULT), which is why it is
        # a separate argument rather than a key of `coherent`.
        cfg.coherent_t_max = float(coherent_t_max)
    if cfg.n_events or cfg.lumi_pb:
        cfg.validate()   # exclusive-and-one-of; leave an unset config alone
    return cfg


def make_plan(name="tensor-thirds", j=1.0, pz=0.7, pzz=0.6, pe=0.7,
              theta_s=0.0, phi_s=None, rel_lumi_offset=0.0, share_plus=0.5,
              pzz_mode="ladder"):
    """A `RunPlan` from a name of `PLANS`.

    tensor-thirds / azz          `tensor_thirds_plan` -- spin 1 only
    helicity-flip / apar         `helicity_flip_plan` at spin `j`
    transverse-tensor / cos2phi  `transverse_tensor_plan` -- spin 1 only
    tensor-flip / flip           `tensor_flip_plan` -- spin 1 only

    `pzz_mode` is a value of `PZZ_MODES` and is READ BY `helicity-flip` /
    `apar` ONLY.  "ladder" (the DEFAULT) is the max-entropy fill at `pz` with
    `pzz` not read at all -- bit for bit what this function did before the
    mode existed; "typed" honours `pzz` through `spin1_populations` /
    `spin32_populations`, refusing a value outside the plan's domain with the
    edge named rather than clamping to it.  THIS IS THE ONE PLACE the mode
    string becomes `HelicityFlipOptions.use_explicit_pzz`.

    THREE OF THE FOUR ARE SPIN-1 PATTERNS AND ARE NOW REFUSED AT ANY OTHER J.
    `tensor_thirds_plan`, `transverse_tensor_plan` and `tensor_flip_plan`
    hard-code `SpinCategory(..., 1.0, spin1_populations(...), ...)` in
    `src/core/bookkeeping.cpp`; only `helicity_flip_plan` takes `j`.  Before
    2026-09-04 the last two BUILT a spin-1 plan at J = 3/2 and the run then
    threw three frames down out of `InclusiveKernel::amplitudes` ("spin state
    J = 1.000000 is not the kernel's ion spin 1.500000"), and
    `tensor-thirds`' own refusal advised `transverse-tensor`, which is one of
    the two that could not work (`phase_D_li7_rank2.md` sec. 1.5, defect F2).
    """
    name = name.lower()
    if pzz_mode not in PZZ_MODES:
        raise ValueError("unknown pzz_mode %r; know %s"
                         % (pzz_mode, ", ".join(PZZ_MODES)))
    if name in ("tensor-thirds", "azz", "transverse-tensor", "cos2phi",
                "tensor-flip", "flip") and abs(j - 1.0) > 1e-12:
        raise ValueError(
            "%s is a SPIN-1 pattern (bookkeeping.cpp hard-codes j = 1 in its "
            "categories) and J = %g was asked for; at J = 3/2 the only plan "
            "that takes j is helicity-flip, and it is a VECTOR plan (the "
            "beam helicity flips, the fill does not).  THERE IS NO SPIN-3/2 "
            "TENSOR PLAN IN THIS TREE, and on the INCLUSIVE channel there "
            "would be nothing for one to measure: 7Li is spin 3/2, "
            "default_inclusive_kernel fills a rank-2 slot for spin 1 only, "
            "so b1_32 = b2_32 = delta_32 = 0 and a T = +1 against T = -1 "
            "contrast comes back as the SAME double (measured at the "
            "default 7Li beam config 1 and scenario: 590952.42641509 pb for "
            "BOTH categories, and the asymmetry exactly 0.0).  Build the two-category T contrast by hand if you "
            "want the zero on the record, use --isotope 6Li for an inclusive "
            "tensor programme, or run --channel tagged-7Li-alpha, whose "
            "alpha-t alignment IS carried (in the event weight, gated at "
            "<P2(cos theta_k)> = -T/5).  See docs/OPEN_ITEMS_SOLUTIONS.md "
            "open item 15." % (name, j))
    if name in ("tensor-thirds", "azz"):
        return _lipolgen.tensor_thirds_plan(pz, pzz, rel_lumi_offset, theta_s,
                                            0.0 if phi_s is None else phi_s)
    if name in ("helicity-flip", "apar"):
        opt = _lipolgen.HelicityFlipOptions()
        opt.theta_s = theta_s
        opt.phi_s = 0.0 if phi_s is None else phi_s
        opt.rel_lumi_offset = rel_lumi_offset
        # THE ONE PLACE `--pzz-mode` becomes a fill.  "ladder" leaves
        # `use_explicit_pzz` false, which is the branch that was always here,
        # so the default is bit for bit by construction and not by an
        # argument about two code paths agreeing.
        if pzz_mode == "typed":
            opt.use_explicit_pzz = True
            opt.pzz = pzz
        return _lipolgen.helicity_flip_plan(j, pz, pe, opt)
    if name in ("transverse-tensor", "cos2phi"):
        return _lipolgen.transverse_tensor_plan(
            pzz, 0.0 if phi_s is None else phi_s)
    if name in ("tensor-flip", "flip"):
        import math
        return _lipolgen.tensor_flip_plan(
            pzz, math.pi / 2.0 if phi_s is None else phi_s, share_plus,
            rel_lumi_offset)
    raise ValueError("unknown plan %r; know %s" % (name, ", ".join(PLANS)))


def make_pipeline(config=None, plan=None, **kw):
    """`Pipeline` from a config/plan or from `make_config` keywords."""
    if config is None:
        config = make_config(**kw)
    if plan is None:
        plan = make_plan("tensor-thirds",
                         j=ion_spin(config.isotope))
    return _lipolgen.Pipeline(config, plan)


def run(isotope="6Li", config=1, channel="inclusive", plan="tensor-thirds",
        events=100000, lumi_pb=0.0, seed=20260713, run_number=1,
        optics="yr-high-acceptance", pz=0.7, pzz=0.6, pe=0.7,
        theta_s=0.0, phi_s=None, pzz_mode="ladder", rel_lumi_offset=0.0,
        nthreads=1, keep_events=False, hadronize=False,
        pythia_options=None, **kw):
    """Build a run and generate it; returns the columnar dict.

    The `pipeline` object is attached to the result under the "pipeline" key
    so that cross sections, optics and the far-forward re-routing stay
    reachable.  `hadronize=True` binds a `PythiaBridge` (T2 tier) into the
    configuration, which forces single-threaded generation.

    THE RUN-PLAN ARGUMENTS ARE THE CLI'S.  `pzz_mode`, `rel_lumi_offset`,
    `theta_s`, `phi_s` and the run number reached neither `make_plan` nor the
    `KnobRunContext` until 2026-09-16: they fell through `**kw` into
    `make_config`, which raised `TypeError: make_config() got an unexpected
    keyword argument 'pzz_mode'`, and `meta["knob_provenance"]` came back with
    68 rows against the CLI's 69 -- the missing one being `pzz_mode` itself.
    The two entry points now write the same table for the same run.
    """
    cfg = make_config(isotope=isotope, config=config, channel=channel,
                      events=events, lumi_pb=lumi_pb, seed=seed,
                      run=run_number, optics=optics, **kw)
    rp = make_plan(plan, j=ion_spin(cfg.isotope), pz=pz, pzz=pzz, pe=pe,
                   theta_s=theta_s, phi_s=phi_s,
                   rel_lumi_offset=rel_lumi_offset, pzz_mode=pzz_mode)
    bridge = None
    if hadronize:
        if not _lipolgen.HAVE_PYTHIA8:
            raise RuntimeError("LiPolGen was built without PYTHIA 8")
        beams = _lipolgen.default_configs(cfg.isotope)[cfg.beam_config]
        opts = _lipolgen.PythiaBridgeOptions()
        for k, v in (pythia_options or {}).items():
            setattr(opts, k, v)
        # ONE unpolarised backend for the whole run (docs/CONVENTIONS.md), the
        # same wiring `cli.py` does: the T2 struck-nucleon SPECIES draw reads
        # `PythiaBridgeOptions.f2_source`, and the `Pipeline` cannot fill it
        # because the bridge is built out here.  Left unset it would keep
        # drawing from `ToyF2` while T0 and T1 moved to the selected backend,
        # and the toy's F2n/F2p is up to 24 % away from CT18NLO's at x = 0.5.
        # None (the `unpol_sf="toy"` default) leaves the bridge's own ToyF2 in
        # place, so this is bit for bit on a default run.  An explicit
        # `pythia_options={"f2_source": ...}` still wins HERE, and then
        # `set_pythia_hadronizer` refuses it unless it is the config's own
        # object: a bridge on a backend `meta["unpol_sf"]` does not name would
        # draw its T2 species from one fit while the rate came from another
        # (up to 24 % on F2n/F2p).  Pass `unpol_sf=` instead of overriding
        # this one field.
        if opts.f2_source is None:
            opts.f2_source = cfg.unpol_sf_obj
        bridge = _lipolgen.PythiaBridge(beams, opts)
        _lipolgen.set_pythia_hadronizer(cfg, bridge)
        nthreads = 1               # PythiaBridge is not re-entrant
    p = _lipolgen.Pipeline(cfg, rp)
    # THE RUN CONTEXT, which this function has always held and never passed.
    # `RunPlan` records its MOMENTS and not which flags produced them, so
    # without a `KnobRunContext` the `pz` / `pzz` / `rel_lumi_offset` rows are
    # omitted rather than guessed and `meta["knob_provenance"]` came back with
    # 63 rows: `lg.run(plan="helicity-flip", pzz=0.6)` recorded the typed
    # `pzz` NOWHERE, which is the silence class the table exists to close, on
    # the library's own entry point.  The T2 half is deliberately left alone:
    # the metadata writer overwrites it from the bridge ACTUALLY bound to the
    # config, which is the only honest source for it.
    ctx = _lipolgen.KnobRunContext()
    ctx.plan_name = plan
    ctx.pz = float(pz)
    ctx.pzz = float(pzz)
    ctx.pe = float(pe)
    ctx.rel_lumi_offset = float(rel_lumi_offset)
    # The row the CLI had and this function did not: without it the table came
    # back with 68 rows against `lipolgen-run`'s 69.
    ctx.pzz_mode = pzz_mode
    # The luminosity is already on the config, so the pipeline resolved the
    # per-category counts in its constructor: generate(0) is the whole run in
    # both modes.
    out = p.generate(0, keep_events, nthreads, context=ctx)
    out["pipeline"] = p
    if bridge is not None:
        out["bridge"] = bridge
    return out
