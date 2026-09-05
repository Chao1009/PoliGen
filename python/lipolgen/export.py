# SPDX-License-Identifier: GPL-3.0-or-later
"""numpy exporters in the schemas `PolarizedLithiumSim/evgen/polligen` uses.

Three schemas, one per consumer:

`inclusive_dict`   the dict `polligen.sample.InclusiveSampler.sample_category`
                   returns -- x, q2, y, phi, m (+ cell where it exists) plus
                   the per-event `category` / `lam_e` labels.

`tagged_dict`      the dict `polligen.tagged.TaggedSampler.sample_category`
                   returns -- the DIS block, the (M, m_S) spin labels, the
                   spectator's spin-frame (k, cos theta_k, phi_k) and its lab
                   kinematics (pT, theta, p_lab, R, xL, kx, ky, kz, phi_spec)
                   and the far-forward `route`.  `polligen.tagged.rp_accepted`
                   reads exactly this dict.

`write_hfs_npz`    the .npz `polligen.hfs.HFSSample.load` reads -- offsets,
                   pid, charge, p4 (E, px, py, pz), x, q2, y, kp, weight,
                   e_energy, p_per_nucleon and a JSON `meta` carrying
                   `sigma_gen_mb`.  The scattered electron is NOT in the
                   particle list (it is `kp`), which is the convention of
                   `tools/pythia8/gen_dis_hfs.py`.

Every function takes either the columnar dict `Pipeline.generate()` returns
(the fast path -- the arrays are already there) or a sequence of `Event`
records (the general path).

`cell` -- the flat accepted-(x, Q2)-cell index `sample.weights_for` needs for
Mode-W reweighting -- is on BOTH paths since 2026-08-30: it is stored on the
record as `Event.kin.cell` (the ION-level sampler's cell on the inclusive
channel, the STRUCK-CLUSTER sampler's on the tagged ones), so a
`Pipeline.generate()` sample reweights exactly like an
`InclusiveSampler.sample_n()` batch.

The tagged dicts also carry the tier-T1 block: `n_partner`, `struck_pdg`,
`struck_pol` and `struck_virtuality`, from the `Role.StruckNucleon` and
`Role.PartnerSpectator` particles the cluster breakup writes (`breakup.hpp`).
These are LiPolGen additions with no polligen counterpart, so `tagged_dict`
puts them in only when they are asked for through `extra` -- except that
`columns_from_events` always builds them.
"""

import json

import numpy as np

from . import _lipolgen as _l

__all__ = ["inclusive_dict", "tagged_dict", "columns_from_events",
           "hfs_sample", "write_hfs_npz", "write_columns_npz",
           "load_hfs_npz", "RC_KEYS", "rc_columns"]

#: `polligen.sample.InclusiveSampler.sample_category` keys.
INCLUSIVE_KEYS = ("x", "q2", "y", "phi", "m", "cell", "category", "lam_e")

#: The tier-T1 columns (`breakup.hpp`); LiPolGen additions, not polligen keys.
T1_KEYS = ("n_partner", "struck_pdg", "struck_pol", "struck_virtuality")

#: rc.hpp's weight columns.  PRESENT ONLY when the run had `--rc` on -- an
#: `--rc off` sample carries exactly today's key set, on both the columnar and
#: the `Event`-record path.  They are NOT in `INCLUSIVE_KEYS` / `TAGGED_KEYS`:
#: those are polligen schemas and do not move.
RC_KEYS = ("rc_tensor_lo", "rc_tensor_hi", "rc_tail")

#: `polligen.tagged.TaggedSampler.sample_category` keys.
TAGGED_KEYS = ("x", "q2", "y", "phi", "m_ion", "m_struck", "k", "cos_theta_k",
               "phi_k", "pT", "theta", "p_lab", "R", "xL", "kx", "ky", "kz",
               "phi_spec", "route", "category", "lam_e")


def _is_columns(obj):
    return isinstance(obj, dict) and "x" in obj


# ------------------------------------------------- Event sequence -> columns

def columns_from_events(events, optics=None, pot_config="18x275"):
    """Columnar dict from a sequence of `Event` records.

    The spectator lab quantities are READ OFF `Event.kin` (`spec_pt`,
    `spec_theta`, `spec_p_lab`, `spec_r`, `spec_xl`, `spec_kx/ky/kz`,
    `phi_spec`) -- `boost_spectator`'s own numbers, stored on the record since
    2026-08-30 instead of being re-derived here.  The re-derivation is exact
    (the beam boost is longitudinal, so (kx, ky) = (px, py) and
    kz = gamma p_z - gamma*beta E inverts the boost) and is kept as a test,
    `python/tests/test_pipeline.py`.  `route` needs an envelope; pass `optics`
    (an `Optics`, e.g. `pipeline.optics`) to get it, otherwise every event is
    `Route.Lost`.
    """
    events = list(events)
    n = len(events)
    out = {k: np.zeros(n) for k in
           ("x", "q2", "y", "phi", "w2", "nu", "s", "weight", "m_ion",
            "m_struck", "pe", "theta_s", "phi_s", "pz", "pzz", "j", "k",
            "cos_theta_k", "phi_k", "alpha_s", "pt_s", "t", "x_pom",
            "e_prime", "theta_e", "eta_e", "pT", "theta", "p_lab", "R", "xL",
            "kx", "ky", "kz", "phi_spec")}
    out["R"][:] = np.nan
    out["xL"][:] = np.nan
    out["struck_virtuality"] = np.full(n, np.nan)
    out["struck_pol"] = np.full(n, 9.0)
    out["kp"] = np.zeros((n, 4))
    out["lam_e"] = np.zeros(n, dtype=np.int64)
    out["route"] = np.zeros(n, dtype=np.int64)
    out["number"] = np.zeros(n, dtype=np.int64)
    out["cell"] = np.full(n, -1, dtype=np.int64)
    out["n_partner"] = np.zeros(n, dtype=np.int64)
    out["struck_pdg"] = np.zeros(n, dtype=np.int64)
    # rc.hpp: the three RC columns exist only when the events carry them, i.e.
    # only when the run had `--rc` on.  The C++ columnar path
    # (`bindings.cpp::columns_to_dict`) applies the SAME rule; missing it here
    # would make `rc_columns` hand back ones for a run that had RC on.
    with_rc = bool(events) and len(events[0].rc_weights) >= len(RC_KEYS)
    if with_rc:
        for key in RC_KEYS:
            out[key] = np.ones(n)
    cats = []
    for i, ev in enumerate(events):
        kin, spin = ev.kin, ev.spin
        for k in ("x", "q2", "y", "phi", "w2", "nu", "s", "k", "cos_theta_k",
                  "phi_k", "alpha_s", "pt_s", "t", "x_pom"):
            out[k][i] = getattr(kin, k)
        out["weight"][i] = ev.weight
        out["cell"][i] = kin.cell
        out["m_ion"][i] = spin.m_ion
        out["m_struck"][i] = spin.m_struck
        for k in ("pe", "theta_s", "phi_s", "pz", "pzz", "j"):
            out[k][i] = getattr(spin, k)
        out["lam_e"][i] = spin.lam_e
        out["number"][i] = ev.number
        if with_rc:
            # Slot 0 -- the event's own PURE spin state.
            rw = ev.rc_weights
            for j, key in enumerate(RC_KEYS):
                out[key][i] = rw[j]
        cats.append(spin.category)
        e1 = ev.find(_l.Role.ScatteredElectron)
        if e1 is not None:
            out["kp"][i] = (e1.p.e, e1.p.px, e1.p.py, e1.p.pz)
            out["e_prime"][i] = e1.p.e
            out["theta_e"][i] = np.arctan2(e1.p.pt(), e1.p.pz)
            pm = e1.p.p()
            if pm > abs(e1.p.pz):
                out["eta_e"][i] = 0.5 * np.log((pm + e1.p.pz) / (pm - e1.p.pz))
        # the tier-T1 block
        pn = ev.find(_l.Role.StruckNucleon)
        if pn is not None:
            out["struck_pdg"][i] = pn.pdg
            out["struck_pol"][i] = pn.pol
            out["struck_virtuality"][i] = pn.p.m2() - _l.M_NUCLEON ** 2
        out["n_partner"][i] = sum(
            1 for q in ev.particles if q.role == _l.Role.PartnerSpectator)

        frag = ev.find(_l.Role.Spectator) or ev.find(_l.Role.IntactRecoil)
        if frag is None:
            continue
        # Stored, not re-derived (see the docstring).
        out["pT"][i] = kin.spec_pt
        out["theta"][i] = kin.spec_theta
        out["p_lab"][i] = kin.spec_p_lab
        out["phi_spec"][i] = kin.phi_spec
        out["kx"][i] = kin.spec_kx
        out["ky"][i] = kin.spec_ky
        out["kz"][i] = kin.spec_kz
        out["R"][i] = kin.spec_r
        out["xL"][i] = kin.spec_xl
        if optics is not None:
            out["route"][i] = _l.route_of(ev, optics, pot_config)
    out["m"] = out["m_ion"]
    out["category"] = np.asarray(cats)
    out["category_names"] = sorted(set(cats))
    return out


def rc_columns(columns):
    """`(rc_tensor_lo, rc_tensor_hi, rc_tail)` arrays from a sample.

    Exactly 1.0 everywhere when the run had `--rc off` (the three keys are then
    absent by design, so an analysis can multiply them in unconditionally).
    Accepts either the columnar dict `Pipeline.generate()` returns or a
    sequence of `Event` records.

    THE WEIGHTS ARE NOT ON `weight`.  `rc_tensor_lo` / `rc_tensor_hi` are a
    two-sided SYSTEMATIC BAND on the tensor part of the rate and `rc_tail` is
    a radiative-tail BACKGROUND; an analysis multiplies them in on purpose:

        w_hi = columns["weight"] * hi        # the +delta edge
        w_bg = columns["weight"] * tail      # Born + radiative tails

    Never quote one band edge alone (`docs/USAGE.md`).
    """
    cols = _columns(columns)
    n = len(cols["x"])
    return tuple(np.asarray(cols[k]) if k in cols else np.ones(n)
                 for k in RC_KEYS)


def _pdg_mass_number(pdg):
    a = abs(int(pdg))
    if a in (2212, 2112):
        return 1
    if a > 1000000000:
        return (a // 10) % 1000
    return 0


def _columns(events, optics=None, pot_config="18x275"):
    if _is_columns(events):
        return events
    return columns_from_events(events, optics, pot_config)


# ----------------------------------------------------------- polligen dicts

def inclusive_dict(events, extra=()):
    """`polligen.sample.InclusiveSampler.sample_category`'s dict.

    Keys x, q2, y, phi, m, cell and the per-event category / lam_e labels.
    Both input paths carry `cell` (`Event.kin.cell` since 2026-08-30), so a
    `Pipeline.generate()` sample reweights like an `InclusiveSampler.sample_n`
    batch.  `extra` names further columns to carry through.
    """
    c = _columns(events)
    out = {}
    for k in INCLUSIVE_KEYS:
        if k in c:
            out[k] = c[k]
    if "m" not in out and "m_ion" in c:
        out["m"] = c["m_ion"]
    for k in extra:
        out[k] = c[k]
    return out


def tagged_dict(events, optics=None, pot_config="18x275", extra=()):
    """`polligen.tagged.TaggedSampler.sample_category`'s dict.

    `polligen.tagged.rp_accepted(events)` consumes this verbatim.  Pass
    `optics` (`pipeline.optics`) when `events` is a sequence of `Event`
    records, so that the far-forward route can be recomputed; a columnar dict
    from `Pipeline.generate()` already carries `route`.
    """
    c = _columns(events, optics, pot_config)
    out = {}
    for k in TAGGED_KEYS:
        if k in c:
            out[k] = c[k]
    if "m_ion" in c:
        out["m"] = c["m_ion"]           # the inclusive spelling, for reuse
    for k in extra:
        out[k] = c[k]
    missing = [k for k in TAGGED_KEYS if k not in out]
    if missing:
        raise KeyError("tagged_dict: input carries none of %s" % missing)
    return out


def rp_accepted(events):
    """`polligen.tagged.rp_accepted`: Roman-Pot main window + near-beam tail."""
    route = np.asarray(events["route"])
    return (route == _l.Route.RomanPots) | (route == _l.Route.RPNearBeam)


# ---------------------------------------------------------------- HFS .npz

def hfs_sample(events, e_energy=None, p_per_nucleon=None, meta=None,
               include_spectators=False, include_hadronic_x=False):
    """The `polligen.hfs.HFSSample` field set, as plain numpy arrays.

    `events` must be a sequence of `Event` records that went through the T2
    (PYTHIA) tier: the particle list is built from `Role.Hadron`,
    `Status.Final` particles.  The scattered electron is excluded and carried
    in `kp`; spectator fragments are excluded too by default, because the
    hadronic-sum identity Sum(E - p_z) = (k + p_N - k')^- that polligen's
    reconstruction is written against is a gamma*-nucleon statement and a
    far-forward alpha would swamp it.
    """
    events = list(events)
    if not events:
        raise ValueError("hfs_sample: no events")
    a = _l.hfs_arrays(events, include_spectators, include_hadronic_x)
    if e_energy is None or p_per_nucleon is None:
        be = events[0].find(_l.Role.BeamElectron)
        bi = events[0].find(_l.Role.BeamIon)
        if be is None or bi is None:
            raise ValueError("hfs_sample: pass e_energy / p_per_nucleon "
                             "(the events carry no beam particles)")
        if e_energy is None:
            e_energy = float(be.p.e)
        if p_per_nucleon is None:
            a_beam = _pdg_mass_number(bi.pdg) or 1
            p_per_nucleon = float(bi.p.pz) / a_beam
    a["e_energy"] = float(e_energy)
    a["p_per_nucleon"] = float(p_per_nucleon)
    a["meta"] = dict(meta or {})
    a["meta"].setdefault("generator", "LiPolGen")
    a["meta"].setdefault("version", _l.__version__)
    a["meta"].setdefault(
        "frame", "head-on: ion +z, electron -z; per-nucleon x,y,Q2")
    return a


def write_hfs_npz(events, path, e_energy=None, p_per_nucleon=None, meta=None,
                  include_spectators=False, include_hadronic_x=False):
    """Write a `polligen.hfs.HFSSample` .npz.

    Read it back with

        from polligen.hfs import HFSSample
        s = HFSSample.load(path)

    `meta` is merged into the JSON blob and should carry `sigma_gen_mb` when
    the sample is meant to be luminosity-merged with another one
    (`HFSSample.concatenate` reweights on it); `Pipeline.generate()["meta"]`
    already has it.
    """
    s = hfs_sample(events, e_energy, p_per_nucleon, meta, include_spectators,
                   include_hadronic_x)
    np.savez_compressed(
        path, offsets=s["offsets"], pid=s["pid"], charge=s["charge"],
        p4=s["p4"], x=s["x"], q2=s["q2"], y=s["y"], kp=s["kp"],
        weight=s["weight"], e_energy=np.array(s["e_energy"]),
        p_per_nucleon=np.array(s["p_per_nucleon"]),
        meta=np.array(json.dumps(_jsonable(s["meta"]))))
    return s


def load_hfs_npz(path):
    """Read a `write_hfs_npz` file back without needing polligen installed."""
    with np.load(path, allow_pickle=False) as f:
        out = {k: f[k] for k in f.files if k != "meta"}
        out["meta"] = json.loads(str(f["meta"])) if "meta" in f.files else {}
    out["e_energy"] = float(out["e_energy"])
    out["p_per_nucleon"] = float(out["p_per_nucleon"])
    return out


# --------------------------------------------------------- columnar .npz

def write_columns_npz(columns, path, keys=None):
    """Write a `Pipeline.generate()` dict to a .npz.

    Non-array entries (`meta`, `category_names`, `pipeline`) are dropped;
    `meta` is written as a JSON string under "meta", the way
    `HFSSample.save` does it.
    """
    arrays = {}
    for k, v in columns.items():
        if k in ("meta", "pipeline", "bridge", "events", "category_names"):
            continue
        v = np.asarray(v)
        if v.dtype.kind in "biufU":
            arrays[k] = v
    if "category_names" in columns:
        arrays["category_names"] = np.asarray(columns["category_names"])
    meta = dict(columns.get("meta", {}))
    arrays["meta"] = np.array(json.dumps(_jsonable(meta)))
    np.savez_compressed(path, **arrays)
    return path


def _jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): _jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    if isinstance(obj, (np.integer,)):
        return int(obj)
    if isinstance(obj, (np.floating,)):
        return float(obj)
    return obj
