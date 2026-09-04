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
    "B1_MODELS", "B1_UNPOL",
    "make_config", "make_plan", "make_pipeline", "run", "__version__",
]

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
#: below x ~ 0.1, COMPARABLE above it -- peak |x b1| 1.67e-4 against Miller's
#: 4.27e-4 at Q2 = 2.5, and 4x LARGER at x = 0.3 with the opposite sign) --
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
#: It reaches the b1 and NOTHING else: the kernel's own `f2_source` stays
#: `ToyF2` on every setting, so the SPIN-BLIND cell cross section
#: (`InclusiveSampler.cell_xsec_pb`) is bit for bit under this flag and only
#: the tensor shift moves.  Read only by "li6-convolution"; `validate()` refuses it
#: on the other two models rather than let `meta["b1_unpol"]` record a PDF
#: that never touched the rate.  `B1UnpolSource.Custom` is not here on
#: purpose: it names an object, so it is set by assigning
#: `config.b1_unpol_sf`.
B1_UNPOL = {
    "toy": _lipolgen.B1UnpolSource.Toy,
    "mstw": _lipolgen.B1UnpolSource.Mstw,
    "ct18nlo": _lipolgen.B1UnpolSource.Ct18Nlo,
}

#: Run-plan names accepted on the command line (aliases included).
PLANS = ("tensor-thirds", "azz", "helicity-flip", "apar", "transverse-tensor",
         "cos2phi", "tensor-flip", "flip")

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
                rc_fq_scale=None, rc_tail_tensor_scale=None,
                rc_qe_suppression=None,
                b1_model=None, b1_band_scale=None,
                b1_alpha_d_dwave_weight=None, b1_unpol=None):
    """A `PipelineConfig` from plain values (the CLI's own constructor).

    `channel` is a key of `CHANNELS`; `optics` a key of `OPTICS`.  The isotope
    a channel implies wins over the `isotope` argument (a 6Li alpha tag is a
    6Li run whatever the caller said), which is `channel_isotope`'s rule.
    `coherent` is a dict of `CoherentScenario` fields.  `cluster_wave` is a
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
    `RcMode` directly, and the five `rc_*` knobs write the matching
    `RcOptions` fields.  The RC weights land on `Event.rc_weights` and on
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
    `b1_unpol` is a key of `B1_UNPOL` ("toy", the default, "mstw" or
    "ct18nlo") or a `B1UnpolSource` directly: the UNPOLARISED backend the
    "li6-convolution" b1 folds its own F1 against.  "mstw" is CDKS's own
    MSTW2008 LO and is the configuration the A = 2 gate passes on; it needs
    the optional PYTHIA tier and raises RuntimeError naming it (or naming the
    missing grid file) in a build that has none -- it is never silently
    replaced by the toy.  It is read only by "li6-convolution" and moves b1
    only: the spin-blind cell cross section is bit-identical across settings,
    the tensor-weighted per-category cross sections are not.
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
                                   rc_fq_scale, rc_tail_tensor_scale,
                                   rc_qe_suppression)):
        opt = cfg.rc_options
        if rc_delta_low_x is not None:
            opt.delta_low_x = float(rc_delta_low_x)
        if rc_delta_high_x is not None:
            opt.delta_high_x = float(rc_delta_high_x)
        if rc_fq_scale is not None:
            opt.fq_scale = float(rc_fq_scale)
        if rc_tail_tensor_scale is not None:
            opt.tail_tensor_scale = float(rc_tail_tensor_scale)
        if rc_qe_suppression is not None:
            opt.qe_suppression = float(rc_qe_suppression)
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
    if cfg.n_events or cfg.lumi_pb:
        cfg.validate()   # exclusive-and-one-of; leave an unset config alone
    return cfg


def make_plan(name="tensor-thirds", j=1.0, pz=0.7, pzz=0.6, pe=0.7,
              theta_s=0.0, phi_s=None, rel_lumi_offset=0.0, share_plus=0.5):
    """A `RunPlan` from a name of `PLANS`.

    tensor-thirds / azz          `tensor_thirds_plan` -- spin 1 only
    helicity-flip / apar         `helicity_flip_plan` at spin `j`
    transverse-tensor / cos2phi  `transverse_tensor_plan`
    tensor-flip / flip           `tensor_flip_plan`
    """
    name = name.lower()
    if name in ("tensor-thirds", "azz"):
        if abs(j - 1.0) > 1e-12:
            raise ValueError(
                "tensor-thirds is the spin-1 equal-thirds pattern; J = %g "
                "needs helicity-flip or transverse-tensor" % j)
        return _lipolgen.tensor_thirds_plan(pz, pzz, rel_lumi_offset, theta_s,
                                            0.0 if phi_s is None else phi_s)
    if name in ("helicity-flip", "apar"):
        opt = _lipolgen.HelicityFlipOptions()
        opt.theta_s = theta_s
        opt.phi_s = 0.0 if phi_s is None else phi_s
        opt.rel_lumi_offset = rel_lumi_offset
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
        events=100000, lumi_pb=0.0, seed=20260713, optics="yr-high-acceptance",
        pz=0.7, pzz=0.6, pe=0.7, nthreads=1, keep_events=False, hadronize=False,
        pythia_options=None, **kw):
    """Build a run and generate it; returns the columnar dict.

    The `pipeline` object is attached to the result under the "pipeline" key
    so that cross sections, optics and the far-forward re-routing stay
    reachable.  `hadronize=True` binds a `PythiaBridge` (T2 tier) into the
    configuration, which forces single-threaded generation.
    """
    cfg = make_config(isotope=isotope, config=config, channel=channel,
                      events=events, lumi_pb=lumi_pb, seed=seed, optics=optics,
                      **kw)
    rp = make_plan(plan, j=ion_spin(cfg.isotope), pz=pz, pzz=pzz, pe=pe)
    bridge = None
    if hadronize:
        if not _lipolgen.HAVE_PYTHIA8:
            raise RuntimeError("LiPolGen was built without PYTHIA 8")
        beams = _lipolgen.default_configs(cfg.isotope)[cfg.beam_config]
        opts = _lipolgen.PythiaBridgeOptions()
        for k, v in (pythia_options or {}).items():
            setattr(opts, k, v)
        bridge = _lipolgen.PythiaBridge(beams, opts)
        _lipolgen.set_pythia_hadronizer(cfg, bridge)
        nthreads = 1               # PythiaBridge is not re-entrant
    p = _lipolgen.Pipeline(cfg, rp)
    # The luminosity is already on the config, so the pipeline resolved the
    # per-category counts in its constructor: generate(0) is the whole run in
    # both modes.
    out = p.generate(0, keep_events, nthreads)
    out["pipeline"] = p
    if bridge is not None:
        out["bridge"] = bridge
    return out
