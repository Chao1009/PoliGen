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

ONE POLLIGEN KEY HAS NO LiPolGen EQUIVALENT ON THE PIPELINE PATH: `cell`, the
flat accepted-(x, Q2)-cell index that `sample.weights_for` needs for Mode-W
reweighting.  It is an `InclusiveSampler` internal and is not stored on the
`Event` record, so it is present in `InclusiveSampler.sample_n()`'s dict and
absent from `Pipeline.generate()`'s.  Reweight sampler batches, not pipeline
output.
"""

import json

import numpy as np

from . import _lipolgen as _l

__all__ = ["inclusive_dict", "tagged_dict", "columns_from_events",
           "hfs_sample", "write_hfs_npz", "write_columns_npz",
           "load_hfs_npz"]

#: `polligen.sample.InclusiveSampler.sample_category` keys.
INCLUSIVE_KEYS = ("x", "q2", "y", "phi", "m", "cell", "category", "lam_e")

#: `polligen.tagged.TaggedSampler.sample_category` keys.
TAGGED_KEYS = ("x", "q2", "y", "phi", "m_ion", "m_struck", "k", "cos_theta_k",
               "phi_k", "pT", "theta", "p_lab", "R", "xL", "kx", "ky", "kz",
               "phi_spec", "route", "category", "lam_e")


def _is_columns(obj):
    return isinstance(obj, dict) and "x" in obj


# ------------------------------------------------- Event sequence -> columns

def columns_from_events(events, optics=None, pot_config="18x275"):
    """Columnar dict from a sequence of `Event` records.

    The spectator lab quantities are RECONSTRUCTED from the record, not read
    off it: the beam boost is longitudinal, so (kx, ky) = (px, py) exactly and
    kz = gamma p_z - gamma*beta E is the exact inverse of `boost_spectator`
    (tagged.hpp).  `route` needs an envelope; pass `optics` (an `Optics`, e.g.
    `pipeline.optics`) to get it, otherwise every event is `Route.Lost`.
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
    out["kp"] = np.zeros((n, 4))
    out["lam_e"] = np.zeros(n, dtype=np.int64)
    out["route"] = np.zeros(n, dtype=np.int64)
    out["number"] = np.zeros(n, dtype=np.int64)
    cats = []
    for i, ev in enumerate(events):
        kin, spin = ev.kin, ev.spin
        for k in ("x", "q2", "y", "phi", "w2", "nu", "s", "k", "cos_theta_k",
                  "phi_k", "alpha_s", "pt_s", "t", "x_pom"):
            out[k][i] = getattr(kin, k)
        out["weight"][i] = ev.weight
        out["m_ion"][i] = spin.m_ion
        out["m_struck"][i] = spin.m_struck
        for k in ("pe", "theta_s", "phi_s", "pz", "pzz", "j"):
            out[k][i] = getattr(spin, k)
        out["lam_e"][i] = spin.lam_e
        out["number"][i] = ev.number
        cats.append(spin.category)
        e1 = ev.find(_l.Role.ScatteredElectron)
        if e1 is not None:
            out["kp"][i] = (e1.p.e, e1.p.px, e1.p.py, e1.p.pz)
            out["e_prime"][i] = e1.p.e
            out["theta_e"][i] = np.arctan2(e1.p.pt(), e1.p.pz)
            pm = e1.p.p()
            if pm > abs(e1.p.pz):
                out["eta_e"][i] = 0.5 * np.log((pm + e1.p.pz) / (pm - e1.p.pz))
        frag = ev.find(_l.Role.Spectator) or ev.find(_l.Role.IntactRecoil)
        if frag is None:
            continue
        bi = ev.find(_l.Role.BeamIon)
        pt, plab = frag.p.pt(), frag.p.p()
        out["pT"][i] = pt
        out["theta"][i] = np.arctan2(pt, frag.p.pz)
        out["p_lab"][i] = plab
        out["phi_spec"][i] = np.arctan2(frag.p.py, frag.p.px)
        out["kx"][i] = frag.p.px
        out["ky"][i] = frag.p.py
        if bi is not None and bi.mass > 0.0:
            gamma, gbeta = bi.p.e / bi.mass, bi.p.pz / bi.mass
            out["kz"][i] = gamma * frag.p.pz - gbeta * frag.p.e
            if frag.charge != 0.0 and bi.charge > 0.0:
                out["R"][i] = (plab / frag.charge) / (bi.p.pz / bi.charge)
            a_frag, a_beam = _pdg_mass_number(frag.pdg), _pdg_mass_number(bi.pdg)
            if a_frag and a_beam:
                out["xL"][i] = plab / (a_frag * bi.p.pz / a_beam)
        if optics is not None:
            out["route"][i] = _l.route_of(ev, optics, pot_config)
    out["m"] = out["m_ion"]
    out["category"] = np.asarray(cats)
    out["category_names"] = sorted(set(cats))
    return out


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

    Keys x, q2, y, phi, m and the per-event category / lam_e labels; `cell`
    only when the input carries it (an `InclusiveSampler.sample_n` batch does,
    a `Pipeline.generate()` sample does not -- see the module docstring).
    `extra` names further columns to carry through.
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
