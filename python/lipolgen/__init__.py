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

from ._lipolgen import *          # noqa: F401,F403
from . import _lipolgen           # noqa: F401
from ._lipolgen import __version__

from . import export              # noqa: F401
from .export import (inclusive_dict, tagged_dict, hfs_sample,  # noqa: F401
                     write_hfs_npz, write_columns_npz)

__all__ = [n for n in dir(_lipolgen) if not n.startswith("_")] + [
    "export", "inclusive_dict", "tagged_dict", "hfs_sample", "write_hfs_npz",
    "write_columns_npz", "CHANNELS", "PLANS", "make_config", "make_plan",
    "make_pipeline", "run", "__version__",
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
                hadronize_coherent=None, apply_optics_lumi_fraction=None,
                cluster_wave=None):
    """A `PipelineConfig` from plain values (the CLI's own constructor).

    `channel` is a key of `CHANNELS`; `optics` a key of `OPTICS`.  The isotope
    a channel implies wins over the `isotope` argument (a 6Li alpha tag is a
    6Li run whatever the caller said), which is `channel_isotope`'s rule.
    `coherent` is a dict of `CoherentScenario` fields.  `cluster_wave` is a
    key of `CLUSTER_WAVES` ("hulthen", the default, or "vmc") or a
    `ClusterWaveSource` directly; "vmc" replaces the lithium alpha-tag radial
    forms with the ANL VMC tables and then ignores `cluster_beta` / `p_d`.
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
    if scenario is not None:
        cfg.scenario = scenario
    if grid is not None:
        cfg.grid = grid
    if hadronize_coherent is not None:
        # The C++ refuses a hadronizer on the coherent channel (C4): PythiaBridge
        # v0 has no coherent-diffractive target and would invent a nucleon that
        # is not in the record's balance.  This is the deliberate escape hatch.
        cfg.hadronize_coherent = bool(hadronize_coherent)
    if apply_optics_lumi_fraction is not None:
        cfg.apply_optics_lumi_fraction = bool(apply_optics_lumi_fraction)
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
