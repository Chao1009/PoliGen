"""THE MATRIX: every knob, on every channel, under every plan -- and the
`knob_provenance` table checked against the OUTPUT HASH of the two runs.

WHY THIS FILE EXISTS.  The repository's rule is that *a knob that did not run
may not be recorded in the `meta` or printed in the banner as if it had*.  It
was enforced knob by knob -- `validate()` refuses `b1_band_scale` on the
Miller branch and `rc_qe_tensor_scale` without the quasi-elastic tail,
`rank2_input_report` labels the 7Li rank-2 zero, `pol_sf_is_read` labelled the
coherent channel -- and it broke FIVE TIMES IN ONE RUN, each time on an axis
the previous fix had not looked at:

  1. the CHANNEL      `--pol-sf` on coherent-6Li (labelled 2026-09-04);
  2. the RUN PLAN     `--pol-sf` under every unpolarised-beam plan, i.e. the
                      CLI's own default `tensor-thirds`, where every category
                      is built at lam_e = 0 and g1 is multiplied by zero --
                      and `--pe` itself, which those plans never read;
  3. the T2 TIER      `--pom-set` / `--pom-rescale` / `--coherent-t2` on the
                      six non-coherent channels, where the Pomeron PYTHIA
                      instance is built and hadronizes ZERO events;
  4. the RC SUB-KNOBS `--rc-fq-scale` and five others on the tagged channels,
                      where `rc_tail_applies` is false BY CONSTRUCTION, and
                      every rc knob on coherent, where `rc_applies` is;
  5. SILENCE          `--cluster-wave`, `--triton-sf`, `--inclusive-b1`,
                      `--cluster-beta`, `--p-d`, `--fsi-sigma-mb` at
                      `--fsi off`, `--coherent-t-max` off the coherent
                      channel -- accepted, some of them MOVING the output, and
                      recorded NOWHERE.

Each round fixed an instance.  This file tests the CLASS, against the one
mechanism that replaced them (`Pipeline::knob_provenance`, pipeline.hpp):

    for every (channel, plan) spec and every knob at a non-default value,
    generate two small runs and assert

      the output hash MOVED       ->  the table says `read`
      the output hash did NOT     ->  the table says `not-read`
      the config was REFUSED      ->  the table says `refused`

    ... and, in the other direction, that EVERY row of the table is covered by
    a variant here or is excused by name with a reason (`test_every_knob_row_
    is_covered`).  A knob added without a matrix entry fails that test, which
    is the whole point: the sixth round cannot be silent.

THE HASH is the one the 2026-09-05 audit used: sha256 over every ndarray
column of `Pipeline.generate` plus `sigma_pb()` and `sigma_per_category_pb()`.
It is deliberately blind to the `meta` -- a test that compared metadata to
metadata would pass on a run where both were wrong.

RUNTIME.  60 events a run and one thread; the baselines are cached per
(spec, base configuration), so the ~40 variants of a spec share four or five
of them.  The whole file is a few tens of seconds.
"""

import hashlib

import numpy as np
import pytest

import lipolgen as lg
from lipolgen import _lipolgen as _l

N = 60
SEED = 7

NO_LHAPDF = ("needs the OPTIONAL LHAPDF tier: --unpol-sf ct18nlo, --pol-sf "
             "nnpdfpol and a caller-supplied kernel on a real fit cannot be "
             "constructed in a build without it (rebuild with "
             "-DLIPOLGEN_WITH_LHAPDF=ON).  The reach RULE for those knobs is "
             "then untested on this build, and the skip says so rather than "
             "passing silently.")
needs_pythia = pytest.mark.skipif(
    not _l.HAVE_PYTHIA8,
    reason="needs the OPTIONAL PYTHIA 8 tier: the T2 knobs --coherent-t2 / "
           "--pom-set / --pom-rescale have no bridge to configure without it "
           "(rebuild with -DLIPOLGEN_WITH_PYTHIA=ON)")


# --------------------------------------------------------------- the specs
#
# `full` specs sweep every variant; the others exist for ONE axis -- the run
# plan -- and sweep only the knobs whose reach can depend on it.  The seven
# full ones are the seven (isotope, channel) combinations the library builds.

class Spec:
    def __init__(self, tag, isotope, channel, plan_name, j,
                 pz=0.7, pzz=0.6, pe=0.7, full=True):
        self.tag, self.isotope, self.channel = tag, isotope, channel
        self.plan_name, self.j = plan_name, j
        self.pz, self.pzz, self.pe, self.full = pz, pzz, pe, full

    def plan(self, pz=None, pzz=None, pe=None, rel=0.0):
        return lg.make_plan(self.plan_name, j=self.j,
                            pz=self.pz if pz is None else pz,
                            pzz=self.pzz if pzz is None else pzz,
                            pe=self.pe if pe is None else pe,
                            rel_lumi_offset=rel)

    def context(self, pz=None, pzz=None, pe=None, rel=0.0):
        c = _l.KnobRunContext()
        c.plan_name = self.plan_name
        c.pz = self.pz if pz is None else pz
        c.pzz = self.pzz if pzz is None else pzz
        c.pe = self.pe if pe is None else pe
        c.rel_lumi_offset = rel
        return c


SPECS = [
    Spec("inclusive-6Li/tensor-thirds", "6Li", "inclusive", "tensor-thirds", 1.0),
    Spec("inclusive-7Li/helicity-flip", "7Li", "inclusive", "helicity-flip", 1.5),
    Spec("inclusive-d/tensor-thirds", "d", "inclusive", "tensor-thirds", 1.0),
    Spec("coherent/tensor-thirds", "6Li", "coherent", "tensor-thirds", 1.0),
    Spec("tagged-6Li-alpha/tensor-thirds", "6Li", "tagged-6Li-alpha",
         "tensor-thirds", 1.0),
    Spec("tagged-7Li-alpha/helicity-flip", "7Li", "tagged-7Li-alpha",
         "helicity-flip", 1.5),
    Spec("tagged-d-p/tensor-thirds", "d", "tagged-d-p", "tensor-thirds", 1.0),
    # THE PLAN AXIS.  `transverse-tensor` and `tensor-flip` take neither --pz
    # nor --pe; `helicity-flip` takes neither --pzz; and `helicity-flip` at
    # --pe 0 carries lam_e != 0 with lam_e * P_e == 0, which is the case that
    # separates the `pe` rule from the `pol_sf` one.
    Spec("inclusive-6Li/helicity-flip", "6Li", "inclusive", "helicity-flip",
         1.0, full=False),
    Spec("inclusive-6Li/transverse-tensor", "6Li", "inclusive",
         "transverse-tensor", 1.0, full=False),
    Spec("inclusive-6Li/tensor-flip", "6Li", "inclusive", "tensor-flip", 1.0,
         full=False),
    Spec("inclusive-6Li/helicity-flip pe=0", "6Li", "inclusive",
         "helicity-flip", 1.0, pe=0.0, full=False),
    Spec("tagged-7Li-alpha/helicity-flip pe=0", "7Li", "tagged-7Li-alpha",
         "helicity-flip", 1.5, pe=0.0, full=False),
]
#: which knobs the plan-axis specs sweep
PLAN_AXIS_KNOBS = ("pz", "pzz", "pe", "rel_lumi_offset", "pol_sf", "unpol_sf",
                   "rc", "rc_delta_low_x", "rc_fq_scale", "rc_scope")


# ------------------------------------------------------------- the variants
#
# One entry per (knob, base configuration): the base kwargs/post-hook, the
# variant kwargs/post-hook, and an optional plan override.  Exactly ONE knob
# differs between base and variant, so a refusal is always attributable to
# the knob the row is named for.

def _post(**attrs):
    def f(cfg):
        for k, v in attrs.items():
            obj, _, field = k.rpartition(".")
            target = cfg
            for part in obj.split(".") if obj else []:
                target = getattr(target, part)
            setattr(target, field, v)
    return f


def _x_max(v):
    def f(cfg):
        sc = cfg.scenario
        sc.x_max = v
        cfg.scenario = sc
    return f


def _chain(*fs):
    def f(cfg):
        for g in fs:
            if g:
                g(cfg)
    return f


def _kernel_ct18(cfg):
    o = _l.InclusiveKernel.Options()
    o.f2_source = _l.LhapdfSF("CT18NLO", 0)
    cfg.kernel = _l.InclusiveKernel(_l.ion_by_name(cfg.isotope), o)


class V:
    """One matrix cell: knob, base, variant.

    `alts` makes the cell ANY-OF instead of one-value: a list of
    (label, kwargs) whose kwargs are each run against the same baseline, and
    the row is expected to be `read` iff AT LEAST ONE of them moves the hash.
    It exists for the three ROUTE knobs and for nothing else.  Their reach is
    a per-event classification -- `Optics::clears` is consulted only for a
    NEAR-BEAM fragment and `over_rigid_route` only for an OVER-RIGID one -- so
    "another value would give another file" is a question about the whole set
    of values a run HAS, not about one of them: at 60 events seed 7,
    `--optics tagging` moves the coherent route column and
    `--optics yr-high-divergence` does not.  `Pipeline::knob_provenance` asks
    exactly this question, by re-routing the run's own sample at exactly these
    values (`Pipeline::RouteReach`, src/core/pipeline.cpp), so the cell and
    the table are the same statement measured twice, independently.
    """

    def __init__(self, knob, label, base=None, var=None, base_post=None,
                 post=None, plan=None, lhapdf=False, alts=None):
        self.knob, self.label = knob, label
        self.base = base or {}
        self.var = var or {}
        self.base_post, self.post = base_post, post
        self.plan = plan or {}
        self.lhapdf = lhapdf
        self.alts = alts

    @property
    def id(self):
        return "%s:%s" % (self.knob, self.label)


XM = _x_max(0.95)
#: The values the ROUTE knobs are swept over, and the ones
#: `Pipeline::RouteReach` re-routes this run's own sample at.  ONE list per
#: knob, here and in src/core/pipeline.cpp; a value in one and not the other
#: would make the cell and the table answer different questions.
OPTICS_ALTS = ("yr-high-divergence", "tagging", "tagging-legacy")
NSIGMA_ALTS = (1.0, 3.0, 30.0)
POT_ALTS = ("5x41", "10x100", "18x275")
#: `--rc tensor-band` with POLRAD Eq. (44)'s quasi-elastic tail switched off.
NO_QE = _post(**{"rc_options.with_qe_tail": False})
CDKS = dict(b1_model="cdks")
CONV = dict(b1_model="li6-convolution")
RC = dict(rc="tensor-band")
VMC = dict(cluster_wave="vmc")
FSI_ON = dict(fsi="glauber-cluster")

VARIANTS = [
    # --- beams, statistics, routing -------------------------------------
    V("seed", "7 -> 8", var=dict(seed=SEED + 1)),
    V("run", "1 -> 2", var=dict(run=2)),
    V("events", "60 -> 40", var=dict(events=40)),
    V("beam_config", "1 -> 2", var=dict(config=2)),
    # --- the spin fill ---------------------------------------------------
    V("pz", "0.7 -> 0.5", plan=dict(pz=0.5)),
    V("pzz", "0.6 -> 0.4", plan=dict(pzz=0.4)),
    V("pe", "-> 0.35", plan=dict(pe=0.35)),
    V("rel_lumi_offset", "0 -> 0.02", plan=dict(rel=0.02)),
    # --- beams, statistics, routing: THE THREE ROUTE KNOBS ----------------
    #
    # No longer excused.  They were kept out of this matrix because "whether a
    # given envelope change flips ANY of 60 events is statistics" -- true of
    # ONE envelope change, and the excuse hid a table that said `read` of all
    # three wherever a route exists while 22 cells of the 2026-09-05 matrix
    # did not move.  The values swept here are exactly the ones the core
    # re-routes its own sample at, so the cell is a second, independent
    # measurement of the same question and not a restatement of the rule.
    V("optics", "yr-high-acceptance -> every other tabulated envelope",
      alts=[("-> " + o, dict(optics=o)) for o in OPTICS_ALTS]),
    V("n_sigma", "10 -> 1 / 3 / 30",
      alts=[("-> %g" % n, dict(n_sigma=n)) for n in NSIGMA_ALTS]),
    V("pot_config", "-> every machine configuration",
      alts=[("-> " + c, dict(pot_config=c)) for c in POT_ALTS]),
    # --- the acceptance window -------------------------------------------
    V("x_max", "1.0 -> 0.95", post=XM),
    # THE TWO COHERENT EDGES, and they are here because the rule was wrong at
    # BOTH of them while 0.95 passed.  0.05 removes 320 of the coherent
    # channel's 1726 rate-carrying cells and moves sigma_pb from 12476.025113
    # to 12446.109030 pb, and the row said NOT READ because the edge was
    # measured on the already-clipped run; 0.10 is bit-identical to the
    # default and the row said READ because exp(log 0.1) = 0.10000000000000002
    # is one ulp above 0.1.  Off the coherent channel both clip the inclusive
    # window and are plainly read, which is the other half of the cell.
    V("x_max", "1.0 -> 0.05", post=_x_max(0.05)),
    V("x_max", "1.0 -> 0.10", post=_x_max(0.10)),
    # inside (largest carrying centre, 0.1]: bit-identical to the default, so
    # the row MUST say not-read -- the cell boundary rule called it read.
    V("x_max", "1.0 -> 0.097", post=_x_max(0.097)),
    V("x_max", "1.0 -> 0.0955", post=_x_max(0.0955)),
    V("x_max", "1.0 -> 0.09549", post=_x_max(0.09549)),
    # --- structure-function selectors ------------------------------------
    V("unpol_sf", "toy -> ct18nlo", var=dict(unpol_sf="ct18nlo"), lhapdf=True),
    V("pol_sf", "toy -> nnpdfpol", var=dict(pol_sf="nnpdfpol"), lhapdf=True),
    V("kernel", "null -> InclusiveKernel(F2 = CT18NLO)", post=_kernel_ct18,
      lhapdf=True),
    V("b1_model", "miller -> cdks (x_max .95)", var=CDKS, base_post=XM,
      post=XM),
    V("b1_band_scale", "1 -> 2 on the default b1_model", var=dict(b1_band_scale=2.0)),
    V("b1_band_scale", "1 -> 2 on cdks (x_max .95)", base=CDKS,
      var=dict(CDKS, b1_band_scale=2.0), base_post=XM, post=XM),
    V("b1_alpha_d_dwave_weight", "1 -> 2 on the default b1_model",
      var=dict(b1_alpha_d_dwave_weight=2.0)),
    V("b1_alpha_d_dwave_weight", "1 -> 2 on li6-convolution (x_max .95)",
      base=CONV, var=dict(CONV, b1_alpha_d_dwave_weight=2.0),
      base_post=XM, post=XM),
    V("b1_unpol", "toy -> ct18nlo on the default b1_model",
      var=dict(b1_unpol="ct18nlo"), lhapdf=True),
    V("b1_unpol", "toy -> ct18nlo on li6-convolution (x_max .95)", base=CONV,
      var=dict(CONV, b1_unpol="ct18nlo"), base_post=XM, post=XM, lhapdf=True),
    # --- the tagged cluster ----------------------------------------------
    V("cluster_wave", "hulthen -> vmc", var=VMC),
    V("cluster_vmc_mc_sigma", "0 -> 1 on hulthen",
      post=_post(cluster_vmc_mc_sigma=1.0)),
    V("cluster_vmc_mc_sigma", "0 -> 1 on vmc", base=VMC, var=VMC,
      post=_post(cluster_vmc_mc_sigma=1.0)),
    V("cluster_beta", "0.30 -> 0.33", var=dict(cluster_beta=0.33)),
    V("cluster_beta", "0.30 -> 0.33 on vmc", base=VMC,
      var=dict(VMC, cluster_beta=0.33)),
    V("p_d", "0.0867 -> 0.10", var=dict(p_d=0.10)),
    V("p_d", "0.0867 -> 0.10 on vmc", base=VMC, var=dict(VMC, p_d=0.10)),
    V("triton_sf", "hulthen -> ciofi-simula", var=dict(triton_sf="ciofi-simula")),
    V("tier", "T1 -> T0", post=_post(tier=_l.Tier.T0)),
    V("inclusive_b1", "off -> on", var=dict(inclusive_b1=True)),
    # --- FSI -------------------------------------------------------------
    V("fsi", "off -> glauber-cluster", var=FSI_ON),
    V("fsi_sigma_mb", "40 -> 20 at fsi off", var=dict(fsi_sigma_mb=20.0)),
    V("fsi_sigma_mb", "40 -> 20 on glauber-cluster", base=FSI_ON,
      var=dict(FSI_ON, fsi_sigma_mb=20.0)),
    # --- RC: the mode, then the two families -----------------------------
    V("rc", "off -> tensor-band", var=RC),
    V("rc_delta_low_x", "0.30 -> 0.19", base=RC,
      var=dict(RC, rc_delta_low_x=0.19)),
    V("rc_delta_high_x", "0.015 -> 0.05", base=RC,
      var=dict(RC, rc_delta_high_x=0.05)),
    V("rc_x_low", "0.01 -> 0.02", base=RC, var=RC,
      post=_post(**{"rc_options.x_low": 0.02})),
    V("rc_x_high", "0.16 -> 0.30", base=RC, var=RC,
      post=_post(**{"rc_options.x_high": 0.30})),
    V("rc_a_transfer_frac", "0 -> 0.5", base=RC,
      var=dict(RC, rc_a_transfer_frac=0.5)),
    V("rc_band_tau_max", "1 -> 1e-3", base=RC, var=RC,
      post=_post(**{"rc_options.band_tau_max": 1e-3})),
    V("rc_scope", "tensor-rate -> tensor-all", base=RC, var=RC,
      post=_post(**{"rc_options.scope": _l.RcScope.TensorAll})),
    V("rc_fq_scale", "1 -> 2", base=RC, var=dict(RC, rc_fq_scale=2.0)),
    V("rc_tail_tensor_scale", "1 -> 2", base=RC,
      var=dict(RC, rc_tail_tensor_scale=2.0)),
    V("rc_c0_shape", "ho -> vmc-ft", base=RC, var=dict(RC, rc_c0_shape="vmc-ft")),
    V("rc_qe_suppression", "1 -> 0.5", base=RC,
      var=dict(RC, rc_qe_suppression=0.5)),
    V("rc_qe_tensor_scale", "0 -> 1", base=RC,
      var=dict(RC, rc_qe_tensor_scale=1.0)),
    V("rc_qe_kf_gev", "0.169 -> 0.25", base=RC, var=RC,
      post=_post(**{"rc_options.qe_kf_gev": 0.25})),
    # THE QUASI-ELASTIC SUB-FAMILY UNDER `with_qe_tail = false`, a base this
    # matrix did not have.  With POLRAD Eq. (44)'s quasi-elastic piece switched
    # off, `qe_suppression` multiplies nothing and `qe_kf_gev` suppresses
    # nothing -- both were accepted, bit-identical, and recorded read /
    # non-default until 2026-09-05, which is exactly what `qe_tensor_scale`
    # has been refused for since the tail knob existed.  All three are refused
    # there now, and these cells are what says so.  (On the tagged and
    # coherent channels the BASE itself is refused -- `with_qe_tail` is a tail
    # knob and the tail does not apply -- so the cell skips and says why.)
    V("rc_qe_suppression", "1 -> 0.5 at with_qe_tail = false", base=RC,
      var=RC, base_post=NO_QE,
      post=_chain(NO_QE, _post(**{"rc_options.qe_suppression": 0.5}))),
    V("rc_qe_kf_gev", "0.169 -> 0.25 at with_qe_tail = false", base=RC,
      var=RC, base_post=NO_QE,
      post=_chain(NO_QE, _post(**{"rc_options.qe_kf_gev": 0.25}))),
    V("rc_qe_tensor_scale", "0 -> 1 at with_qe_tail = false", base=RC,
      var=RC, base_post=NO_QE,
      post=_chain(NO_QE, _post(**{"rc_options.qe_tensor_scale": 1.0}))),
    V("rc_tail_model", "t-peak -> t-peak+ll", base=RC,
      var=dict(RC, rc_tail_model="t-peak+ll")),
    V("rc_with_qe_tail", "true -> false", base=RC, var=RC,
      post=_post(**{"rc_options.with_qe_tail": False})),
    V("rc_with_tail", "true -> false", base=RC, var=RC,
      post=_post(**{"rc_options.with_tail": False})),
    V("rc_n_eta", "128 -> 64", base=RC, var=RC,
      post=_post(**{"rc_options.n_eta": 64})),
    # 0.01, not 0.5: the ceiling is on the tail RATIO and the shipped tail is
    # only a few per cent, so a 0.5 ceiling clips nothing at 60 events (it
    # does at 200).  The test is about whether the knob is CONSULTED.
    V("rc_tail_max", "10 -> 0.01 (clips every event)", base=RC, var=RC,
      post=_post(**{"rc_options.tail_max": 0.01})),
    V("rc_m_lepton", "M_ELECTRON -> the muon mass", base=RC, var=RC,
      post=_post(**{"rc_options.m_lepton": 0.1056583755})),
    # --- the coherent channel --------------------------------------------
    V("coherent_t_max", "0.2 -> 0.1", var=dict(coherent_t_max=0.1)),
    V("coherent_f0", "0.04 -> 0.05", var=dict(coherent=dict(f0=0.05))),
    V("coherent_slope_b", "50 -> 40", var=dict(coherent=dict(slope_b=40.0))),
    V("coherent_amp", "0.01 -> 0.02", var=dict(coherent=dict(amp=0.02))),
    V("coherent_eps_b0", "-0.08 -> -0.05",
      var=dict(coherent=dict(eps_b0=-0.05))),
    V("coherent_m_x_min", "1.2 -> 1.5",
      post=_post(**{"coherent_xpom.m_x_min": 1.5})),
    V("coherent_x_pom_max", "0.1 -> 0.05",
      post=_post(**{"coherent_xpom.x_pom_max": 0.05})),
    V("coherent_weighted_azimuth", "false -> true",
      post=_post(coherent_weighted_azimuth=True)),
]

#: Rows the generation matrix cannot vary, each with the reason.  They are
#: NOT skipped silently: `test_every_knob_row_is_covered` requires every row
#: to be either here or in `VARIANTS`, so a new knob has to be classified.
EXCUSED = {
    "isotope": "it IS the spec: varying it is a different run, not a knob "
               "variation of this one",
    "channel": "it IS the spec, and it is what every other row is scoped by",
    "lumi_pb": "reachable only by leaving fixed-count mode, which moves "
               "`events` at the same time -- two knobs, not one",
    "poisson": "luminosity mode only, and reaching it moves `events` too",
    "apply_optics_lumi_fraction":
        "luminosity mode only (it multiplies COUNTS and never a cross "
        "section), and reaching it moves `events` too",
    "with_virtual_photon":
        "it adds a status-3 documentation particle to the EVENT RECORD and "
        "no generated COLUMN, so the column hash cannot see it -- the record "
        "path is tests/test_pipeline.cpp's",
    "hadronize": "a T2 knob: covered by test_the_t2_knobs_are_read_on_the "
                 "_coherent_channel_alone, which needs the PYTHIA tier",
    "coherent_t2": "a T2 knob -- see `hadronize`",
    "pom_set": "a T2 knob -- see `hadronize`",
    "pom_rescale": "a T2 knob -- see `hadronize`",
}


# ------------------------------------------------------------------ driver

def _hash(p, cols):
    h = hashlib.sha256()
    for k in sorted(cols):
        v = cols[k]
        if isinstance(v, np.ndarray):
            h.update(k.encode())
            h.update(str(v.dtype).encode())
            h.update(np.ascontiguousarray(v).tobytes())
    h.update(np.array([p.sigma_pb()]).tobytes())
    h.update(np.array(p.sigma_per_category_pb()).tobytes())
    return h.hexdigest()[:24]


def _run(spec, kw, post, plan_kw):
    """(status, hash, table) of one cell.

    `status` is 'ran' or 'refused'; a refusal from `PipelineConfig::validate`
    and one from the `Pipeline` constructor are the same thing to this test --
    both are the run surface saying no before an event exists.
    """
    args = dict(isotope=spec.isotope, channel=spec.channel, config=1,
                events=N, seed=SEED)
    args.update(kw)
    try:
        cfg = lg.make_config(**args)
        if post:
            post(cfg)
        cfg.validate()
        plan = spec.plan(**plan_kw)
        p = _l.Pipeline(cfg, plan)
    except (RuntimeError, ValueError) as e:
        return "refused", None, None, str(e)
    cols = p.generate(0, False, 1)
    rows = {r.name: r for r in p.knob_provenance(spec.context(**plan_kw))}
    return "ran", _hash(p, cols), rows, ""


@pytest.fixture(scope="module")
def cache():
    return {}


def _baseline(cache, spec, v):
    key = (spec.tag, repr(sorted(v.base.items())), id(v.base_post))
    if key not in cache:
        cache[key] = _run(spec, v.base, v.base_post, {})
    return cache[key]


CELLS = [(s, v) for s in SPECS for v in VARIANTS
         if s.full or v.knob in PLAN_AXIS_KNOBS]


@pytest.mark.parametrize(
    "spec,variant", CELLS, ids=["%s|%s" % (s.tag, v.id) for s, v in CELLS])
def test_the_table_matches_the_output_hash(cache, spec, variant):
    """THE gate.  Moved -> `read`; did not move -> `not-read`; refused ->
    `refused`.  Nothing else is accepted, in either direction."""
    if variant.lhapdf and not _l.HAVE_LHAPDF:
        pytest.skip(NO_LHAPDF)
    base = _baseline(cache, spec, variant)
    if base[0] != "ran":
        pytest.skip("the BASE of this cell is refused on %s, so there is no "
                    "run to vary: %s" % (spec.tag, base[3][:200]))
    rows = base[2]
    assert variant.knob in rows, (
        "%s has no knob_provenance row -- every knob a run accepts must have "
        "one" % variant.knob)
    row = rows[variant.knob]

    if variant.alts:
        # THE ANY-OF CELL (the three route knobs).  `read` promises that SOME
        # other value of this knob gives another file, so the cell has to try
        # the values this run HAS -- a value the run surface refuses (the
        # tagging working point is tabulated for 6Li and 7Li only, so
        # `--optics tagging` throws on a deuteron beam) is not one of them and
        # is skipped, not counted as "did not move".
        moved, first, ran = False, None, []
        for lab, kw in variant.alts:
            got, h, _, msg = _run(spec, dict(variant.var, **kw), variant.post,
                                  variant.plan)
            if got == "refused":
                continue
            ran.append(lab)
            if h != base[1]:
                moved, first = True, lab
                break
        assert ran, ("%s on %s: every alternative value is refused by the run "
                     "surface, so this cell measures nothing"
                     % (variant.knob, spec.tag))
        if moved:
            assert row.status == _l.KnobStatus.Read, (
                "%s on %s: %s MOVED the output hash and the table says %s -- "
                "%s" % (variant.knob, spec.tag, first,
                        _l.knob_status_name(row.status), row.reason[:200]))
        else:
            assert row.status == _l.KnobStatus.NotRead, (
                "%s on %s: every value this run has (%s) is BIT-IDENTICAL to "
                "the baseline in all array columns and every sigma, and the "
                "table says `%s` -- %s"
                % (variant.knob, spec.tag, ", ".join(ran),
                   _l.knob_status_name(row.status), row.reason[:200]))
        return

    got, h, vrows, msg = _run(spec, variant.var, variant.post, variant.plan)
    if got == "refused":
        # A NOT-READ row PROMISES that the value is accepted and labelled --
        # that is what makes a label the honest answer instead of a refusal.
        # If the run surface says no, the row is wrong about the surface.
        # `read` is allowed here: a knob can be consulted at its own value and
        # still have other values refused (b1_model = miller on an inclusive
        # deuteron run reaches the rate; cdks is refused because the transfer
        # is 6Li's).
        assert row.status != _l.KnobStatus.NotRead, (
            "%s on %s: the run surface REFUSES a non-default value (%s) while "
            "the table calls it `not-read`, which promises the value is "
            "accepted and labelled.  Record the axis as `refused`."
            % (variant.knob, spec.tag, msg[:200]))
        return

    # IT IS THE VARIANT'S OWN ROW THAT IS JUDGED FROM HERE ON, not the
    # baseline's.  The two can honestly differ, and the coherent `x_max` cell
    # is why: at the default 1.0 the knob clips nothing and the row correctly
    # says not-read ("a value BELOW the edge would move this channel"), while
    # at 0.05 it removes 320 rate-carrying cells and the row must say read.
    # The measurement here -- "this VALUE moved the file" -- is a statement
    # about the variant, and the `meta` that would carry a misleading value is
    # the variant's.  (The refused branch above keeps the baseline's row: a
    # refused run writes no row at all.)
    row = vrows[variant.knob]
    # It ran, so the axis is NOT closed and `refused` is a false claim.
    assert row.status != _l.KnobStatus.Refused, (
        "%s on %s: the table calls the axis `refused` -- validate() throws on "
        "any other value -- but %s was ACCEPTED and built."
        % (variant.knob, spec.tag, variant.label))

    moved = (h != base[1])
    if moved:
        assert row.status == _l.KnobStatus.Read, (
            "%s on %s: %s MOVED the output hash, and the table says %s -- "
            "%s.  A knob that ran must be recorded as read."
            % (variant.knob, spec.tag, variant.label,
               _l.knob_status_name(row.status), row.reason[:200]))
    else:
        assert row.status == _l.KnobStatus.NotRead, (
            "%s on %s: %s is BIT-IDENTICAL to the baseline in all array "
            "columns and every sigma, and the table says `%s` -- %s.  This "
            "is the defect the table exists to prevent: a knob that did not "
            "run, recorded as if it had."
            % (variant.knob, spec.tag, variant.label,
               _l.knob_status_name(row.status), row.reason[:200]))


# ------------------------------------------------------- the other direction

@pytest.mark.parametrize("spec", SPECS, ids=[s.tag for s in SPECS])
def test_every_knob_row_is_covered(spec):
    """A knob added without a matrix entry FAILS HERE.

    This is what makes the file a class fix rather than five instance fixes:
    the table must enumerate every knob, and every row of it must be either
    exercised above or excused by name with a reason.
    """
    got, _, rows, msg = _run(spec, {}, None, {})
    assert got == "ran", msg
    covered = {v.knob for v in VARIANTS} | set(EXCUSED)
    missing = sorted(set(rows) - covered)
    assert not missing, (
        "knob_provenance rows with no matrix entry and no excuse: %s.  Add a "
        "V(...) to VARIANTS, or an entry to EXCUSED saying why the generation "
        "matrix cannot vary it." % ", ".join(missing))
    stale = sorted(covered - set(rows))
    assert not stale, (
        "matrix entries for knobs the table no longer has: %s" % ", ".join(stale))


@pytest.mark.parametrize("spec", SPECS, ids=[s.tag for s in SPECS])
def test_every_row_is_well_formed(spec):
    """Shape, not physics: a row with no reason is a row that explains
    nothing, and a NOT-READ row with no label is one the `meta` cannot use."""
    got, _, rows, msg = _run(spec, {}, None, {})
    assert got == "ran", msg
    assert len(rows) > 40
    for name, r in rows.items():
        assert r.name == name
        assert r.value != "", name
        assert r.reason.strip() != "", name
        if r.status == _l.KnobStatus.NotRead:
            assert r.label.startswith("not read"), (name, r.label)
            assert r.reason.startswith(r.label), (name, r.label, r.reason)
            assert r.meta_value == r.label, name
        else:
            assert r.meta_value == r.value, name


def test_the_meta_carries_the_whole_table():
    """`meta["knob_provenance"]` is the table, not a summary of it."""
    spec = SPECS[0]
    cfg = lg.make_config(isotope=spec.isotope, channel=spec.channel,
                         events=N, seed=SEED)
    p = _l.Pipeline(cfg, spec.plan())
    ctx = spec.context()
    meta = p.generate(1, False, 1, context=ctx)["meta"]
    rows = p.knob_provenance(ctx)
    prov = meta["knob_provenance"]
    assert set(prov) == {r.name for r in rows}
    for r in rows:
        e = prov[r.name]
        assert e["value"] == r.value
        assert e["status"] == _l.knob_status_name(r.status)
        assert e["reason"] == r.reason
        assert e["at_default"] == r.at_default
    # ... and it survives the npz round trip, which is JSON.
    import json
    from lipolgen import export
    assert json.loads(json.dumps(export._jsonable(dict(meta))))[
        "knob_provenance"]["pol_sf"]["status"] == "not-read"


def test_the_banner_block_is_printed_from_the_same_table():
    from lipolgen import cli
    spec = SPECS[0]
    cfg = lg.make_config(isotope=spec.isotope, channel=spec.channel,
                         events=N, seed=SEED, pol_sf="toy")
    rows = _l.Pipeline(cfg, spec.plan()).knob_provenance(spec.context())
    text = "\n".join(cli.knob_provenance_lines(rows))
    assert "KNOB PROVENANCE" in text
    # every NOT-READ knob is named, and its scope with it
    for r in rows:
        if r.status == _l.KnobStatus.NotRead:
            assert r.name in text, r.name
            assert r.label in text, r.label
    # the reason is spelled out for the ones that would mislead
    loud = [r for r in rows
            if r.status == _l.KnobStatus.NotRead and not r.at_default]
    for r in loud:
        assert r.reason[:60] in text, r.name


# ----------------------------------------------------------- the route knobs

@pytest.mark.parametrize(
    "channel,isotope,plan_name,j",
    [("inclusive", "6Li", "tensor-thirds", 1.0),
     ("inclusive", "d", "tensor-thirds", 1.0),
     ("coherent", "6Li", "tensor-thirds", 1.0),
     ("tagged-6Li-alpha", "6Li", "tensor-thirds", 1.0),
     ("tagged-7Li-alpha", "7Li", "helicity-flip", 1.5),
     ("tagged-d-p", "d", "tensor-thirds", 1.0)])
def test_the_route_knobs_are_read_where_a_route_moves(channel, isotope,
                                                      plan_name, j):
    """The same rule as the matrix above, at 2000 events instead of 60 -- and
    that is the point of running it twice.

    `--optics`, `n_sigma` and `pot_config` reach exactly one quantity, the
    per-event `route` column, and they reach it through exactly two branches
    of `route_charged`: the envelope is consulted only for a NEAR-BEAM
    fragment (|R - 1| < 0.05) and `pot_config` only for an OVER-RIGID one
    (R > 1.05).  `route_of` returns `Route.Lost` BEFORE it looks at either
    when the event carries neither a tagged spectator nor an intact recoil --
    which is every INCLUSIVE event.

    WHAT THE OLD FORM OF THIS TEST GOT WRONG.  It asserted `read` if and only
    if the channel writes a far-forward fragment, which is the CLASSIFIER
    being invoked and not the knob being read: at 2000 events seed 11 NONE of
    the eight alternative envelopes and machine configurations moves a single
    route label on tagged-d-p, and `--pot-config` moves none on the coherent
    channel at any statistics.  The rule is the matrix's rule -- moved ->
    `read`, did not move -> `not-read` -- and it is measured here at a second
    statistics precisely because the answer is allowed to differ between them
    (it does: on tagged-7Li-alpha `--optics` moves nothing at 60 events and
    does move at 2000, and the table follows the run rather than a table of
    channels).
    """
    plan = lg.make_plan(plan_name, j=j, pz=0.7, pzz=0.6, pe=0.7)

    def run(**kw):
        cfg = lg.make_config(isotope=isotope, channel=channel, events=2000,
                             seed=11, **kw)
        p = _l.Pipeline(cfg, plan)
        cols = p.generate(0, False, 1)
        return (np.asarray(cols["route"]),
                {r.name: r for r in p.knob_provenance()})

    base, rows = run()
    has_route = channel != "inclusive"
    alts = {"optics": [dict(optics=o) for o in OPTICS_ALTS],
            "n_sigma": [dict(n_sigma=n) for n in NSIGMA_ALTS],
            "pot_config": [dict(pot_config=c) for c in POT_ALTS]}
    for knob, kws in alts.items():
        moved = False
        for kw in kws:
            try:
                other, _ = run(**kw)
            except RuntimeError:
                continue          # not a value this run has
            if not np.array_equal(base, other):
                moved = True
                break
        assert (rows[knob].status == _l.KnobStatus.Read) == moved, (
            knob, channel, moved, rows[knob].reason[:300])
        # ... and the row says which it is on THIS run, not on channels in
        # general: every reason names the measurement.
        assert "re-rout" in rows[knob].reason or not has_route, knob

    if not has_route:
        # THE MEASUREMENT THE TABLE'S REASON QUOTES, and it is absolute: not
        # "no event happened to move" but "the envelope is never reached".
        assert np.all(base == int(_l.Route.Lost))
        for k in ("optics", "n_sigma", "pot_config"):
            assert rows[k].status == _l.KnobStatus.NotRead
            assert "route` column is 0" in rows[k].reason
    else:
        # ... and on the other side a far-forward fragment IS written and
        # priced, whether or not any envelope moves its label.
        roles = {pa.role for ev in _l.Pipeline(
            lg.make_config(isotope=isotope, channel=channel, events=50,
                           seed=11), plan).generate(0, True, 1)["events"]
            for pa in ev.particles}
        assert roles & {_l.Role.Spectator, _l.Role.IntactRecoil}, channel


def test_pot_config_is_never_read_on_the_coherent_channel():
    """The one route cell that is structural and not statistical.

    `pot_config` enters `route_charged` through `over_rigid_route` alone,
    which is tested at R > 1 + NEAR_BEAM_BAND; the coherent channel's
    far-forward fragment is the INTACT beam nucleus, which can only have LOST
    momentum, so it is never over-rigid and the branch is never reached.  The
    table says so on both statistics and names the branch.
    """
    plan = lg.tensor_thirds_plan(0.7, 0.6)
    for n_events in (60, 2000):
        cfg = lg.make_config(isotope="6Li", channel="coherent",
                             events=n_events, seed=11)
        p = _l.Pipeline(cfg, plan)
        cols = p.generate(0, False, 1)
        r = np.asarray(cols["R"])
        assert np.nanmax(r) <= 1.0 + 1e-9, np.nanmax(r)
        row = {q.name: q for q in p.knob_provenance()}["pot_config"]
        assert row.status == _l.KnobStatus.NotRead
        assert "over_rigid_route" in row.reason
        # ... and the npz key carries the label, not the machine name.
        meta = p.generate(1, False, 1)["meta"]
        assert meta["pot_config"] == row.label != p.pot_config


# ---------------------------------------------------------------- the T2 tier

@needs_pythia
@pytest.mark.parametrize("channel,isotope",
                         [("coherent", "6Li"), ("inclusive", "6Li")])
def test_the_t2_knobs_are_read_on_the_coherent_channel_alone(channel, isotope):
    """Class 3.  The Pomeron PYTHIA instance is built whenever
    `coherent_t2 == Pomeron` REGARDLESS of channel and hadronizes ZERO events
    off the coherent one, while `meta` recorded `pom_set` as if it had run."""
    plan = lg.tensor_thirds_plan(0.7, 0.6)

    def run(pom_set, t2="pomeron", hadronize=True):
        cfg = lg.make_config(isotope=isotope, channel=channel, events=20,
                             seed=4242)
        bridge = None
        if hadronize:
            popts = _l.PythiaBridgeOptions()
            popts.coherent_t2 = (_l.CoherentT2.Pomeron if t2 == "pomeron"
                                 else _l.CoherentT2.Off)
            popts.pom_set = pom_set
            beams = _l.default_configs(cfg.isotope)[cfg.beam_config]
            bridge = _l.PythiaBridge(beams, popts)
            _l.set_pythia_hadronizer(cfg, bridge)
        p = _l.Pipeline(cfg, plan)
        cols = p.generate(0, True, 1)
        pids = [pa.pdg for ev in cols["events"] for pa in ev.particles]
        # THE TABLE THE METADATA WRITER BUILT.  It fills the T2 half of the
        # KnobRunContext from the bridge ACTUALLY bound to the config
        # (`PythiaHadronizerHook`), which is the only honest source for it --
        # a caller-supplied context could say anything.
        meta = dict(cols["meta"])
        rows = meta["knob_provenance"]
        return pids, rows, meta

    a_pid, a_rows, a_meta = run(6)
    b_pid, b_rows, _ = run(5)
    want_read = (channel == "coherent")
    for k in ("pom_set", "pom_rescale", "coherent_t2"):
        assert (a_rows[k]["status"] == "read") == want_read, (
            k, channel, a_rows[k]["reason"][:200])
    if want_read:
        assert a_pid != b_pid          # the set moves the hadrons here
        assert a_meta["pom_set"] == 6  # ... and is recorded as the int it is
    else:
        # MEASURED: not one particle moves, so the label goes in the meta and
        # not the number.
        assert a_pid == b_pid
        assert a_meta["pom_set"] == a_rows["pom_set"]["label"]
        assert a_meta["pom_set"] != 6
        assert a_meta["coherent_t2"] == a_rows["coherent_t2"]["label"]

    # ... and without --hadronize the three are recorded, labelled, rather
    # than nowhere at all (the fourth defect class).
    _, off_rows, off_meta = run(5, hadronize=False)
    for k in ("pom_set", "pom_rescale", "coherent_t2"):
        assert off_rows[k]["status"] == "not-read"
        assert "--hadronize" in off_rows[k]["reason"]
        assert k not in off_meta          # the block stays conditional ...
