#!/usr/bin/env python3
# SPDX-License-Identifier: GPL-3.0-or-later
"""Dump numeric reference tables from the Python `polligen` generator.

This script imports the EXISTING Python event generator `polligen`
(PolarizedLithiumSim/evgen/polligen) -- it does not reimplement any
physics -- and writes plain-JSON tables of the numbers a C++
reimplementation (LiPolGen) must reproduce at rtol 1e-12.

It touches nothing under PolarizedLithiumSim; it only imports it.

Only the TOY backends are exercised (ToyF2, ToyG1, toy_b1(mode="toy"),
toy_delta_gluon, r_sigma_lt): every one of those is a closed-form
formula (no PDF grids, no digitized CSV tables), so a from-scratch C++
port needs only the formulas recorded in each JSON's "toy_formulas"
block plus the grid/argument values recorded alongside the outputs.

Run: python3 dump_polligen_reference.py
Output: LiPolGen/validation/reference/*.json (this script also writes
LiPolGen/validation/README.md, describing every table).
"""

import json
import math
import pathlib
import sys

import numpy as np

HERE = pathlib.Path(__file__).resolve().parent
EVGEN = HERE.parents[1] / "PolarizedLithiumSim" / "evgen"
FASTSIM = HERE.parents[1] / "PolarizedLithiumSim" / "fastsim"
OUT = HERE / "reference"
OUT.mkdir(parents=True, exist_ok=True)

if not EVGEN.is_dir():
    raise SystemExit("evgen package not found at %s" % EVGEN)

sys.path.insert(0, str(EVGEN))
# polligen/__init__.py inserts fastsim/ itself, but do it here too so
# that `from polli_fastsim import ...` works even before `import polligen`.
if str(FASTSIM) not in sys.path:
    sys.path.insert(0, str(FASTSIM))

import polligen  # noqa: E402
from polligen import spin, xsec, bookkeeping, tagged, coherent  # noqa: E402

from polli_fastsim import beams, structure, polarized, spectator  # noqa: E402
from polli_fastsim import asymmetries as asym  # noqa: E402

MISSING = []  # (function, reason) pairs for functions we could not call


# ---------------------------------------------------------------------------
# JSON-safe conversion: everything becomes plain python floats/ints/lists so
# that json.dump's default float formatter (repr(), i.e. shortest
# round-tripping decimal = full double precision) is used verbatim.
# ---------------------------------------------------------------------------

def j(x):
    """Recursively convert numpy / python objects to JSON-safe types."""
    if isinstance(x, dict):
        return {str(k): j(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [j(v) for v in x]
    if isinstance(x, np.ndarray):
        return j(x.tolist())
    if isinstance(x, complex):
        return {"re": float(x.real), "im": float(x.imag)}
    if isinstance(x, (np.complexfloating,)):
        return {"re": float(x.real), "im": float(x.imag)}
    if isinstance(x, (np.floating,)):
        return float(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (np.bool_, bool)):
        return bool(x)
    if isinstance(x, float):
        return x
    if isinstance(x, int):
        return x
    if x is None:
        return None
    if isinstance(x, str):
        return x
    raise TypeError("cannot JSON-ify %r (%s)" % (x, type(x)))


def cplx_matrix(m):
    """A complex numpy matrix -> {"re": [[...]], "im": [[...]]}."""
    m = np.asarray(m)
    return {"re": j(m.real), "im": j(m.imag)}


def write_json(name, obj):
    path = OUT / name
    with open(path, "w") as fh:
        json.dump(j(obj), fh, indent=1, sort_keys=True)
        fh.write("\n")
    return path


README_SECTIONS = []


def doc(title, text):
    README_SECTIONS.append((title, text))


# ---------------------------------------------------------------------------
# Shared constants
# ---------------------------------------------------------------------------

AXES = [
    (0.0, 0.0),
    (np.pi / 2.0, 0.0),
    (np.pi / 2.0, np.pi / 3.0),
    (np.arccos(1.0 / np.sqrt(3.0)), 0.7),
    (1.1, 2.2),
]

XQ2_POINTS = [(0.0224, 1.14), (0.056, 1.14), (0.141, 3.13), (0.141, 14.3),
              (0.005, 2.0), (0.30, 30.0)]

PE = 0.7


# ===========================================================================
# 1. spin.json
# ===========================================================================

def build_spin():
    out = {"axes": [{"theta": t, "phi": p} for t, p in AXES],
           "wigner_d": [], "clebsch_gordan": [], "spins": {}}

    # --- Wigner-d at a few arguments (module: polligen.spin.wigner_d) ------
    for jj, beta in [(1.0, 0.7), (1.5, 1.1), (1.0, 0.3), (1.5, 2.0)]:
        out["wigner_d"].append({
            "j": jj, "beta": beta,
            "m_order": spin.m_values(jj).tolist(),
            "d": spin.wigner_d(jj, beta).tolist(),
        })

    # --- Clebsch-Gordan at a few arguments (polligen.spin.clebsch_gordan) --
    cg_args = [
        (1, 1, 1, -1, 0, 0), (0.5, 0.5, 0.5, -0.5, 1, 0),
        (1, 0, 1, 0, 2, 0), (2, 0, 1, 1, 1, 1), (2, 0, 1, 0, 1, 0),
        (1.5, 1.5, 0.0, 0.0, 1.5, 1.5), (1, 1, 0.5, -0.5, 1.5, 0.5),
        (1, -1, 0.5, 0.5, 0.5, -0.5),
    ]
    for a in cg_args:
        out["clebsch_gordan"].append({"args": list(a),
                                      "value": spin.clebsch_gordan(*a)})

    for jj, label in [(1.0, "j1"), (1.5, "j3_2")]:
        ms = spin.m_values(jj).tolist()
        n = len(ms)
        pops_list = []

        def add(name, pops):
            pops_list.append((name, tuple(float(p) for p in pops)))

        for i, m in enumerate(ms):
            pure = [0.0] * n
            pure[i] = 1.0
            add("pure_m=%g" % m, pure)
        add("uniform", [1.0 / n] * n)

        pz_anchor = 8.0 / 13.0 if abs(jj - 1.0) < 1e-9 else 0.7
        pz_cross = 0.7 if abs(jj - 1.0) < 1e-9 else 8.0 / 13.0
        add("maxent_pz=%.10g" % pz_anchor, spin.populations_maxent(jj, pz_anchor))
        add("maxent_pz=%.10g" % pz_cross, spin.populations_maxent(jj, pz_cross))

        if abs(jj - 1.0) < 1e-9:
            add("explicit_pz=0.7_pzz=0.4", spin.spin1_populations(0.7, 0.4))
            add("explicit_pz=8/13_pzz=4/13",
                spin.spin1_populations(8.0 / 13.0, 4.0 / 13.0))
        else:
            add("explicit_pz=0.7_t=0.4", spin.spin32_populations(0.7, 0.4, 0.0))
            add("explicit_pz=8/13_t=4/13",
                spin.spin32_populations(8.0 / 13.0, 4.0 / 13.0, 0.0))

        pop_entries = []
        for name, pops in pops_list:
            moments_own = spin.moments_along_axis(jj, pops)
            entry = {
                "name": name, "populations": list(pops),
                "moments_own_axis": {
                    "vector": moments_own[0], "tensor_zz": moments_own[1],
                } if abs(jj - 1.0) < 1e-9 else {
                    "vector": moments_own[0], "tensor_zz": moments_own[1],
                    "octupole_z": moments_own[2],
                },
                "axes": [],
            }
            for theta, phi in AXES:
                dens = spin.SpinDensity(jj, pops, theta, phi)
                rho = dens.rho_lab()
                mom = dens.lab_moments()
                axis_entry = {
                    "theta": theta, "phi": phi,
                    "rho": cplx_matrix(rho),
                    "vector_polarization": spin.vector_polarization(rho, jj).tolist(),
                    "tensor_polarization": spin.tensor_polarization(rho, jj),
                    "lab_moments": {
                        "vector": mom["vector"].tolist(),
                        "tensor_zz": mom["tensor_zz"],
                    },
                }
                if abs(jj - 1.5) < 1e-9:
                    axis_entry["octupole_moment"] = spin.octupole_moment(rho, jj)
                    axis_entry["lab_moments"]["octupole_z"] = mom["octupole_z"]
                entry["axes"].append(axis_entry)
            pop_entries.append(entry)

        out["spins"][label] = {"j": jj, "m_order": ms,
                                "population_sets": pop_entries}
    return out


doc("spin.json",
    "Spin-density matrices and multipole moments, from `polligen.spin`\n"
    "(PolarizedLithiumSim/evgen/polligen/spin.py).\n\n"
    "- `wigner_d`: `spin.wigner_d(j, beta)` (spin.py:68), package m-ordering\n"
    "  (+J...-J, see `spin.m_values`, spin.py:40).\n"
    "- `clebsch_gordan`: `spin.clebsch_gordan(j1,m1,j2,m2,j,m)` (spin.py:46),\n"
    "  Racah formula, exact for these half-integer/integer spins.\n"
    "- `spins.j1` / `spins.j3_2`: for J=1 and J=3/2, several population\n"
    "  vectors p_m (m ordered +J...-J):\n"
    "    * `pure_m=<m>`      -- one entry 1, rest 0\n"
    "    * `uniform`         -- all populations equal\n"
    "    * `maxent_pz=<pz>`  -- `spin.populations_maxent(j, pz)` (spin.py:210),\n"
    "      spin-temperature p_m ~ exp(beta m); pz=8/13 (J=1) and pz=0.7 (J=3/2)\n"
    "      are the t=3 rational anchors of evgen/tests/test_bookkeeping.py,\n"
    "      giving tensor moments 4/13 and 0.4 respectively (bisection tol\n"
    "      1e-13, so results are good to ~1e-13, not full double precision --\n"
    "      compare at a looser tolerance than rtol=1e-12 for this table only).\n"
    "    * `explicit_pz=..._pzz=...` / `..._t=...` -- `spin.spin1_populations`\n"
    "      (spin.py:180) / `spin.spin32_populations` (spin.py:193), closed-form,\n"
    "      exact to double precision.\n"
    "  For each population vector: `moments_own_axis` is\n"
    "  `spin.moments_along_axis(j, populations)` (spin.py:165) -- (vector,\n"
    "  tensor_zz[, octupole_z]) along the population's OWN axis. Then for\n"
    "  each of the 5 `axes` (theta, phi): `rho` is\n"
    "  `spin.SpinDensity(j, populations, theta, phi).rho_lab()` (spin.py:247),\n"
    "  a (2J+1)x(2J+1) complex matrix (m-ordering +J...-J), split into\n"
    "  `re`/`im`; `vector_polarization` is `spin.vector_polarization(rho, j)`\n"
    "  (spin.py:138); `tensor_polarization` is `spin.tensor_polarization(rho, j)`\n"
    "  (spin.py:144, J=1: <3Jz^2-2>; J=3/2: <3Jz^2-J(J+1)>/3);\n"
    "  `octupole_moment` (J=3/2 only) is `spin.octupole_moment(rho, j)`\n"
    "  (spin.py:156); `lab_moments` is `SpinDensity.lab_moments()`\n"
    "  (spin.py:251), which must equal (moment_own_axis * n_hat) for the\n"
    "  vector and (tensor_own_axis * P2(cos theta)) for tensor_zz -- this is\n"
    "  the analytic identity `test_rotated_moments_analytic` in\n"
    "  evgen/tests/test_spin.py pins.")


# ===========================================================================
# 2. xsec.json
# ===========================================================================

def kernel_amplitude_dump(kern, ion, s, points, axes, ms, with_target_mass=False,
                          with_tensor_gamma=False):
    xs = np.array([p[0] for p in points])
    q2s = np.array([p[1] for p in points])
    ys = q2s / (s * xs)
    t = kern.tables(xs, q2s, with_g2=True)
    tables_out = {k: v.tolist() for k, v in t.items()}

    r_vals = structure.r_sigma_lt(xs, q2s)
    asym_out = {
        "y": ys.tolist(),
        "R": r_vals.tolist(),
        "depolarization_d": asym.depolarization_d(ys, xs, q2s).tolist(),
        "a_parallel": asym.a_parallel(t["g1"], t["f1"], ys, xs, q2s).tolist(),
    }
    if "b1" in t:
        asym_out["azz"] = asym.azz(t["b1"], t["f1"], t["f2"], xs, ys,
                                    b2=t["b2"]).tolist()
    if "delta" in t:
        asym_out["a_cos2phi"] = asym.a_cos2phi(t["delta"], t["f1"], t["f2"],
                                               xs, ys).tolist()

    amp_entries = []
    for lam_e in (+1, -1):
        for theta_s, phi_s in axes:
            for m in ms:
                state = xsec.EventSpinState(lam_e, PE, ion.spin, m,
                                            theta_s=theta_s, phi_s=phi_s)
                w_avg, a1, a2 = kern.amplitudes(t, xs, q2s, s, state,
                                                with_perp=True)
                amp_entries.append({
                    "lam_e": lam_e, "pe": PE, "m": m,
                    "theta_s": theta_s, "phi_s": phi_s,
                    "w_avg": w_avg.tolist(), "a1": a1.tolist(),
                    "a2": a2.tolist(),
                })
    out = {"tables": tables_out, "asymmetries": asym_out,
           "amplitudes": amp_entries}

    # dsigma at a few phi values for two representative states
    phis = np.linspace(0.0, 2.0 * np.pi, 8, endpoint=False)
    dsig_entries = []
    for lam_e, theta_s, phi_s, m in (
            (+1, 0.0, 0.0, ms[0]),
            (+1, np.pi / 2.0, 0.0, ms[0]),
            (0, np.pi / 2.0, 0.3, ms[0])):
        state = xsec.EventSpinState(lam_e, PE, ion.spin, m, theta_s=theta_s,
                                    phi_s=phi_s)
        rows = []
        for x0, q20 in points:
            dsig = kern.dsigma(np.full(phis.shape, x0),
                               np.full(phis.shape, q20), phis, s, state,
                               with_perp=True)
            rows.append(dsig.tolist())
        dsig_entries.append({
            "lam_e": lam_e, "pe": PE, "m": m, "theta_s": theta_s,
            "phi_s": phi_s, "phi_grid": phis.tolist(),
            "dsigma_per_point": rows,
        })
    out["dsigma"] = dsig_entries
    out["dsigma_unpol"] = kern.dsigma_unpol(xs, q2s, s).tolist()

    if with_target_mass:
        gamma2 = xsec.gamma_squared(xs, q2s)
        eps_g = xsec.epsilon_gamma(ys, gamma2)
        d_gamma = xsec.depolarization_gamma(ys, gamma2, r_vals)
        eta_g = xsec.eta_gamma(ys, gamma2)
        out["target_mass"] = {
            "gamma_squared": gamma2.tolist(),
            "epsilon_gamma": eps_g.tolist(),
            "depolarization_gamma": d_gamma.tolist(),
            "eta_gamma": eta_g.tolist(),
            "a_parallel_finite_gamma": kern.a_parallel(t, xs, q2s, ys).tolist(),
            "a_parallel_massless": asym.a_parallel(t["g1"], t["f1"], ys, xs,
                                                    q2s).tolist(),
        }

    if with_tensor_gamma:
        # The EXACT finite-gamma tensor sector (Cosyn Eqs. 9/10/14/16/17/24,
        # plans/08 D2).  The kernel-level amplitudes above already carry it
        # (this kernel was built with tensor_gamma=True); what is dumped here
        # in addition are the three pieces it is assembled from, so that a
        # port can be located rather than only compared:
        #   theta_q          Eq. (24), the rest-frame photon-beam angle
        #   cosyn_tensor_sfs Eqs. (17a)-(17e) on the block's own b1..b4
        #   cosyn_unpol_sfs  Eq. (16)
        #   harmonics        (h0, h1, h2) per (axis, m)
        gamma2 = xsec.gamma_squared(xs, q2s)
        cq, sq = xsec.theta_q_cos_sin(ys, gamma2)
        f_t, f_l, f_lt, f_tt = xsec.cosyn_tensor_sfs(
            t["b1"], t["b2"], t["b3"], t["b4"], xs, gamma2)
        fu_t, fu_l = xsec.cosyn_unpolarized_sfs(t["f1"], t["f2"], xs, gamma2)
        harm = []
        for theta_s, phi_s in axes:
            for m in ms:
                state = xsec.EventSpinState(0, 0.0, ion.spin, m,
                                            theta_s=theta_s, phi_s=phi_s)
                h0, h1, h2 = kern._tensor_harmonics_gamma(t, xs, q2s, ys,
                                                          state)
                harm.append({
                    "m": m, "theta_s": theta_s, "phi_s": phi_s,
                    "h0": np.broadcast_to(h0, xs.shape).tolist(),
                    "h1": np.broadcast_to(h1, xs.shape).tolist(),
                    "h2": np.broadcast_to(h2, xs.shape).tolist(),
                })
        out["tensor_gamma"] = {
            "cos_theta_q": np.broadcast_to(cq, xs.shape).tolist(),
            "sin_theta_q": np.broadcast_to(sq, xs.shape).tolist(),
            "F_TLL_T": np.broadcast_to(f_t, xs.shape).tolist(),
            "F_TLL_L": np.broadcast_to(f_l, xs.shape).tolist(),
            "F_TLT": np.broadcast_to(f_lt, xs.shape).tolist(),
            "F_TTT": np.broadcast_to(f_tt, xs.shape).tolist(),
            "F_UU_T": np.broadcast_to(fu_t, xs.shape).tolist(),
            "F_UU_L": np.broadcast_to(fu_l, xs.shape).tolist(),
            "harmonics": harm,
        }
    return out


# The higher-twist tensor slots the `kernel_tensor_gamma` blocks are built
# with.  b3 and b4 are UNMEASURED; these are scenario shapes chosen only so
# that both slots reach the kernel and neither cancels the other.
B3_SCENARIO = "b3_func = 0.05 * f1"
B4_SCENARIO = "b4_func = -0.02 * f1"


def build_xsec():
    toy_b1_formula = ("toy_b1(x, q2, f1, mode='toy') = "
                       "0.01 * max(x,1e-6)^-0.2 * (1 - x/0.20) * exp(-3*x) * f1")
    toy_delta_formula = ("toy_delta_gluon(x, q2, f1, scale=1e-3) = "
                         "scale * f1 * x^0.3 * (1-x)^4")
    r_formula = "r_sigma_lt(x, q2) = 0.18 / (1 + q2/50.0)"

    out = {
        "constants": {
            "ALPHA_EM": structure.ALPHA_EM,
            "GEV2_TO_PB": structure.GEV2_TO_PB,
            "TENSOR_LL_SIGN": asym.TENSOR_LL_SIGN,
            "LI6_CLUSTER_POLARIZATION": beams.LI6_CLUSTER_POLARIZATION,
            "M_NUCLEON": xsec.M_NUCLEON,
            "ions": {
                ion.name: {
                    "A": ion.A, "Z": ion.Z, "N": ion.N, "spin": ion.spin,
                    "eff_pol_p": ion.eff_pol_p, "eff_pol_n": ion.eff_pol_n,
                    "mass_per_nucleon": ion.mass_per_nucleon,
                } for ion in (beams.LI6, beams.LI7)
            },
        },
        "toy_formulas": {
            "ToyF2.f2p(x,q2)": ("lq=ln(max(q2,1.1)/0.04); lam=0.045*lq; "
                                "sea=0.20*x^(-lam)*(1-x)^7; "
                                "val=1.05*x^0.55*(1-x)^3; f2p=sea+val"),
            "ToyF2.f2n_over_f2p(x)": "clip(1 - 0.75*x, 0.25, 1.0)",
            "ToyF2.f2n(x,q2)": "f2p(x,q2) * f2n_over_f2p(x)",
            "NuclearF2.f2a(x,q2)": "Z*f2p(x,q2) + N*f2n(x,q2)  (emc_ratio=None)",
            "NuclearF2.f1a(x,q2)": "f2a(x,q2) / (2*x*(1+R)),  R=r_sigma_lt",
            "r_sigma_lt(x,q2)": r_formula,
            "ToyG1.a1p(x)": "clip(x^0.7, 0, 1)",
            "ToyG1.a1n(x)": "-0.07*(1-x)^2 + 0.8*x^2.2",
            "ToyG1.g1p(x,q2)": "a1p(x) * f2p(x,q2)/(2*x*(1+R))",
            "ToyG1.g1n(x,q2)": "a1n(x) * f2n(x,q2)/(2*x*(1+R))",
            "ToyG1.g1_nucleus(ion,x,q2)":
                "Z*eff_pol_p*g1p(x,q2) + N*eff_pol_n*g1n(x,q2)",
            "toy_b1": toy_b1_formula,
            "toy_delta_gluon": toy_delta_formula,
            "b2_default": "b2 = 2*x*b1  (InclusiveKernel.tables, no b2_func)",
            "g2_ww": ("Wandzura-Wilczek: g2(x) = -g1(x) + "
                     "int_x^1 g1(u)/u du (fixed Q2); computed by "
                     "xsec.g2_ww via a u=x^(1-t) substitution, trapezoid, "
                     "npts=96 default"),
        },
        "grid_points": [{"x": x, "q2": q2} for x, q2 in XQ2_POINTS],
        "axes": [{"theta": t, "phi": p} for t, p in AXES],
        "results": {},
    }

    for ion_name, ion, ms in (
            ("6Li", beams.LI6, spin.m_values(1.0).tolist()),
            ("7Li", beams.LI7, spin.m_values(1.5).tolist())):
        config = beams.default_configs(ion_name)[1]
        s = config.sqrt_s_per_nucleon ** 2

        b1_func = lambda x, q2, f1: polarized.toy_b1(x, q2, f1, mode="toy")  # noqa: E731
        delta_func = lambda x, q2, f1: polarized.toy_delta_gluon(x, q2, f1, scale=1e-3)  # noqa: E731
        b3_func = lambda x, q2, f1: 0.05 * f1  # noqa: E731
        b4_func = lambda x, q2, f1: -0.02 * f1  # noqa: E731

        entry = {
            "beam_config": {"electron_energy": config.electron_energy,
                            "ion_momentum_per_nucleon":
                                config.ion_momentum_per_nucleon,
                            "sqrt_s_per_nucleon": config.sqrt_s_per_nucleon,
                            "s": s},
        }

        if abs(ion.spin - 1.0) < 1e-9:
            kern = xsec.InclusiveKernel(ion, b1_func=b1_func,
                                        delta_func=delta_func,
                                        target_mass=False)
            entry["kernel_default"] = kernel_amplitude_dump(
                kern, ion, s, XQ2_POINTS, AXES, ms, with_target_mass=False)
            kern_tm = xsec.InclusiveKernel(ion, b1_func=b1_func,
                                           delta_func=delta_func,
                                           target_mass=True)
            entry["kernel_default_target_mass"] = kernel_amplitude_dump(
                kern_tm, ion, s, XQ2_POINTS, AXES, ms, with_target_mass=True)
            # the EXACT finite-gamma tensor sector, off by default
            kern_tg = xsec.InclusiveKernel(ion, b1_func=b1_func,
                                           delta_func=delta_func,
                                           b3_func=b3_func, b4_func=b4_func,
                                           tensor_gamma=True,
                                           target_mass=False)
            entry["kernel_tensor_gamma"] = kernel_amplitude_dump(
                kern_tg, ion, s, XQ2_POINTS, AXES, ms,
                with_tensor_gamma=True)
            entry["kernel_tensor_gamma"]["toy_formulas"] = {
                "b3_func": B3_SCENARIO, "b4_func": B4_SCENARIO,
            }
        else:
            # rank-0/1 default kernel: rank-2 (b1_32/delta_32) default to
            # zero, exactly as InclusiveKernel documents for spin 3/2.
            kern = xsec.InclusiveKernel(ion, target_mass=False)
            entry["kernel_default"] = kernel_amplitude_dump(
                kern, ion, s, XQ2_POINTS, AXES, ms, with_target_mass=False)
            kern_tm = xsec.InclusiveKernel(ion, target_mass=True)
            entry["kernel_default_target_mass"] = kernel_amplitude_dump(
                kern_tm, ion, s, XQ2_POINTS, AXES, ms, with_target_mass=True)
            # scenario rank-2 kernel (matches
            # evgen/tests/test_tensor_convention.py's spin-3/2 lambdas)
            b1_32 = lambda x, q2, f1: 0.05 * f1  # noqa: E731
            delta_32 = lambda x, q2, f1: -1e-2 * f1  # noqa: E731
            kern2 = xsec.InclusiveKernel(ion, b1_32_func=b1_32,
                                         delta_32_func=delta_32,
                                         target_mass=False)
            entry["kernel_scenario_rank2"] = kernel_amplitude_dump(
                kern2, ion, s, XQ2_POINTS, AXES, ms, with_target_mass=False)
            entry["kernel_scenario_rank2"]["toy_formulas"] = {
                "b1_32_func": "0.05 * f1",
                "delta_32_func": "-1e-2 * f1",
            }
            # the finite-gamma tensor sector on the spin-3/2 rank-2 scenario,
            # which also exercises the b3/b4 slots on a J = 3/2 alignment
            kern2g = xsec.InclusiveKernel(ion, b1_32_func=b1_32,
                                          delta_32_func=delta_32,
                                          b3_func=b3_func, b4_func=b4_func,
                                          tensor_gamma=True,
                                          target_mass=False)
            entry["kernel_tensor_gamma"] = kernel_amplitude_dump(
                kern2g, ion, s, XQ2_POINTS, AXES, ms,
                with_tensor_gamma=True)
            entry["kernel_tensor_gamma"]["toy_formulas"] = {
                "b1_32_func": "0.05 * f1",
                "delta_32_func": "-1e-2 * f1",
                "b3_func": B3_SCENARIO, "b4_func": B4_SCENARIO,
            }
        out["results"][ion_name] = entry
    return out


doc("xsec.json",
    "`polligen.xsec.InclusiveKernel` (PolarizedLithiumSim/evgen/polligen/xsec.py)\n"
    "and the underlying `polli_fastsim.asymmetries` analytic formulas\n"
    "(PolarizedLithiumSim/fastsim/polli_fastsim/asymmetries.py), evaluated\n"
    "on TOY structure functions only.\n\n"
    "`constants`: `structure.ALPHA_EM` (structure.py:21, =1/137.036),\n"
    "`structure.GEV2_TO_PB` (structure.py:38), `asymmetries.TENSOR_LL_SIGN`\n"
    "(=-1.0 since 2026-08-29: the LITERATURE convention, Cosyn et al. Eq. 27\n"
    "/ HERMES, A_zz = -(2/3) b1/F1), `beams.LI6_CLUSTER_POLARIZATION`\n"
    "(=0.81123, the whole-nucleus 6Li vector polarization of the cluster\n"
    "picture, from which beams.LI6's per-nucleon slots are a third each),\n"
    "`asymmetries.M_NUCLEON` (=0.9383 GeV,\n"
    "the FREE-nucleon mass used in gamma^2 = 4 M^2 x^2/Q^2), and per-ion\n"
    "`beams.Ion` fields (A, Z, N=A-Z, spin, eff_pol_p, eff_pol_n,\n"
    "mass_per_nucleon) for `beams.LI6`/`beams.LI7` (beams.py:197,199).\n\n"
    "`toy_formulas`: the exact closed forms used everywhere below --\n"
    "`ToyF2` (structure.py:41-59), `r_sigma_lt` (structure.py:62),\n"
    "`NuclearF2.f2a`/`f1a` (structure.py:303-338), `ToyG1`\n"
    "(polarized.py:26-82), `toy_b1(..., mode='toy')` (polarized.py:444,451 --\n"
    "NOTE: the default mode='digitized' reads a CSV table and is NOT used\n"
    "here; every b1 in this file is the pure formula), `toy_delta_gluon`\n"
    "(polarized.py:536), the b2 default `2*x*b1` (xsec.py `tables()`,\n"
    "xsec.py:234/238), and `g2_ww` (xsec.py:85, Wandzura-Wilczek by\n"
    "quadrature -- reproduce it with the SAME npts=96 trapezoid rule to\n"
    "match at 1e-12; it is not a closed form).\n\n"
    "`results.<6Li|7Li>`: `beam_config` is `beams.default_configs(name)[1]`\n"
    "(beams.py:224, the MID energy point) with its `sqrt_s_per_nucleon`\n"
    "(beams.py:213) and `s = sqrt_s_per_nucleon^2`.\n\n"
    "`kernel_default` (6Li: J=1, b1_func/delta_func = the toy formulas\n"
    "above, target_mass=False; 7Li: J=3/2, all rank-2 kwargs None so\n"
    "b1/b2/delta are identically zero per xsec.py's documented default) and\n"
    "`kernel_default_target_mass` (same kernel with target_mass=True, the\n"
    "exact finite-gamma E143 vector sector) each contain:\n"
    "  * `tables`: `InclusiveKernel.tables(x, q2, with_g2=True)` (xsec.py:217)\n"
    "    -- f1, f2, g1, g2, b1, b2, delta on `grid_points`.\n"
    "  * `asymmetries`: `y = q2/(s x)`, `R = r_sigma_lt(x,q2)`,\n"
    "    `asymmetries.depolarization_d` (asymmetries.py:44),\n"
    "    `asymmetries.a_parallel` (asymmetries.py:59, the MASSLESS g1/F1\n"
    "    form -- this is what `kernel_default`'s vector sector must equal\n"
    "    exactly, per `test_vector_sector_matches_a_parallel`),\n"
    "    `asymmetries.azz` (asymmetries.py:101, only when b1 is nonzero)\n"
    "    and `asymmetries.a_cos2phi` (asymmetries.py:111, only when delta\n"
    "    is nonzero).\n"
    "  * `amplitudes`: `InclusiveKernel.amplitudes(tables, x, q2, s, state,\n"
    "    with_perp=True)` (xsec.py:345) -- (w_avg, a1, a2) of\n"
    "    W = 1 + w_avg + a1 cos(phi') + a2 cos(2 phi') -- for every\n"
    "    (lam_e in {+1,-1}) x (5 `axes`) x (each spin projection m of the\n"
    "    ion), at pe=0.7 (`PE` in this script).\n"
    "  * `dsigma`: `InclusiveKernel.dsigma(x, q2, phi, s, state,\n"
    "    with_perp=True)` (xsec.py:388) on `phi_grid` (8 uniform points in\n"
    "    [0, 2pi)) for 3 representative states, at each of the 6\n"
    "    `grid_points` (`dsigma_per_point[i]` is the phi-array at\n"
    "    grid_points[i]).  `dsigma_unpol` is `InclusiveKernel.dsigma_unpol`\n"
    "    (xsec.py:383), i.e. `structure.dsigma_dx_dq2` (structure.py:341) on\n"
    "    the per-nucleon F2.\n"
    "  * (tensor_gamma variant only) `tensor_gamma`: the pieces of the EXACT\n"
    "    finite-gamma tensor sector (Cosyn et al. Eqs. 9/10/14/16/17/24,\n"
    "    plans/08 D2) -- `xsec.theta_q_cos_sin` (Eq. 24),\n"
    "    `xsec.cosyn_tensor_sfs` (Eqs. 17a-17e, on the block's own b1..b4),\n"
    "    `xsec.cosyn_unpolarized_sfs` (Eq. 16), and the (h0, h1, h2)\n"
    "    harmonics `InclusiveKernel._tensor_harmonics_gamma` returns for\n"
    "    every (axis, m).  Those blocks are dumped from a kernel built with\n"
    "    `tensor_gamma=True` and both higher-twist slots filled\n"
    "    (b3_func=0.05*f1, b4_func=-0.02*f1 -- SCENARIO shapes, b3 and b4\n"
    "    are unmeasured), so their `amplitudes` carry the exact b-sector\n"
    "    while every other block carries the massless one.\n"
    "  * (target_mass variant only) `target_mass`: `xsec.gamma_squared`\n"
    "    (xsec.py:118), `xsec.epsilon_gamma` (xsec.py:129),\n"
    "    `xsec.depolarization_gamma` (xsec.py:136), `xsec.eta_gamma`\n"
    "    (xsec.py:142), `InclusiveKernel.a_parallel` with target_mass=True\n"
    "    (xsec.py:260, the exact E143 D_gamma(A1+eta*A2) form) alongside the\n"
    "    massless `asymmetries.a_parallel` for comparison.\n\n"
    "`results.7Li.kernel_scenario_rank2`: the SAME structure with an\n"
    "explicit rank-2 scenario (b1_32_func=0.05*f1, delta_32_func=-1e-2*f1,\n"
    "matching `evgen/tests/test_tensor_convention.py`'s\n"
    "`test_spin32_rate_and_cos2phi_channels_are_now_consistent`), so the\n"
    "J=3/2 tensor geometry Q_NN=(3m^2-J(J+1))/3, t_geo=Q_NN*P2(cos theta_S),\n"
    "c_eff=3*Q_NN (xsec.py:311-341) is exercised numerically too.")


# ===========================================================================
# 3. beams.json
# ===========================================================================

def build_beams():
    out = {"constants": {
        "PROTON_TOP_MOMENTUM": beams.PROTON_TOP_MOMENTUM,
        "PROTON_MASS": beams.PROTON_MASS,
        "PROTON_CONFIG_ENERGIES": list(beams.PROTON_CONFIG_ENERGIES),
        "ELECTRON_ENERGIES": list(beams.ELECTRON_ENERGIES),
        "NUCLEUS_MASS": {"%s_%d_%d" % k: v
                         for k, v in beams.NUCLEUS_MASS.items()},
        # the 6Li cluster wave function -- ONE source of truth for the
        # inclusive effective polarization below and for the tagged S/D
        # interference of `polligen.tagged`, which re-exports these two
        # (beams.py, author decision 2026-08-29 / plans/04 #6)
        "P_D_LI6": beams.P_D_LI6,
        "P_D_DEUTERON": beams.P_D_DEUTERON,
        "ALPHA_D_VECTOR_POLARIZATION": beams.ALPHA_D_VECTOR_POLARIZATION,
        "DEUTERON_VECTOR_POLARIZATION": beams.DEUTERON_VECTOR_POLARIZATION,
        "LI6_CLUSTER_POLARIZATION": beams.LI6_CLUSTER_POLARIZATION,
        "LI6_NAIVE_ONE_THIRD": beams.LI6_NAIVE_ONE_THIRD,
    }, "ions": {}, "configs": {}}

    for name, ion in beams.IONS.items():
        out["ions"][name] = {
            "A": ion.A, "Z": ion.Z, "N": ion.N, "spin": ion.spin,
            "eff_pol_p": ion.eff_pol_p, "eff_pol_n": ion.eff_pol_n,
            "mass_per_nucleon": ion.mass_per_nucleon,
            "momentum_per_nucleon_max": ion.momentum_per_nucleon_max,
        }

    for name in ("6Li", "7Li", "d"):
        configs = beams.default_configs(name)
        rows = []
        for cfg in configs:
            gamma_ion = math.sqrt(cfg.ion_momentum_per_nucleon ** 2
                                  + cfg.ion.mass_per_nucleon ** 2) \
                / cfg.ion.mass_per_nucleon
            rows.append({
                "electron_energy": cfg.electron_energy,
                "ion_momentum_per_nucleon": cfg.ion_momentum_per_nucleon,
                "sqrt_s_per_nucleon": cfg.sqrt_s_per_nucleon,
                "gamma_ion": gamma_ion,
                "label": cfg.label(),
            })
        out["configs"][name] = rows

    out["gamma_of_proton_energy"] = [
        {"proton_energy": pe, "gamma": beams.gamma_of(pe)}
        for pe in beams.PROTON_CONFIG_ENERGIES
    ]
    return out


doc("beams.json",
    "`polli_fastsim.beams` (PolarizedLithiumSim/fastsim/polli_fastsim/beams.py).\n\n"
    "`constants`: `PROTON_TOP_MOMENTUM` (beams.py:35, =275.0 GeV),\n"
    "`PROTON_MASS` (beams.py:36), `PROTON_CONFIG_ENERGIES` (beams.py:52,\n"
    "the 41/100/275 GeV proton energies), `ELECTRON_ENERGIES` (beams.py:203,\n"
    "5/10/18 GeV), `NUCLEUS_MASS` (beams.py:43, keyed `\"<name>_<A>_<Z>\"`).\n\n"
    "`ions`: every `beams.IONS` entry (beams.py:194-201) -- A, Z, N=A-Z,\n"
    "spin, eff_pol_p, eff_pol_n, `Ion.mass_per_nucleon` (beams.py:169, the\n"
    "PHYSICAL nuclear mass / A, not amu), `Ion.momentum_per_nucleon_max`\n"
    "(beams.py:176, = 275*Z/A GeV).\n\n"
    "`configs.<6Li|7Li|d>`: `beams.default_configs(name)` (beams.py:224),\n"
    "the low/mid/top (5x41-like / 10x100-like / 18x275-like) reference scan,\n"
    "with `BeamConfig.sqrt_s_per_nucleon` (beams.py:213,\n"
    "= sqrt(4*E_e*p_ion_per_nucleon)) and `gamma_ion = "
    "sqrt(p^2+m^2)/m` computed here from `ion_momentum_per_nucleon` and\n"
    "`Ion.mass_per_nucleon` (not stored directly by beams.py, but a plain\n"
    "function of the two fields above).\n\n"
    "`gamma_of_proton_energy`: `beams.gamma_of(proton_energy)` (beams.py:86)\n"
    "at the three reference proton energies -- the RING Lorentz factor that\n"
    "every ion at that configuration is gamma-matched to.")


# ===========================================================================
# 4. tagged.json
# ===========================================================================

def wave_dump(channel_name, channel):
    k = np.linspace(1e-4, 1.2, 40)
    waves = []
    for w in channel.waves:
        waves.append({
            "l": w.l, "prob": w.prob, "beta": w.beta,
            "kappa": channel.base.kappa,
            "k_grid": k.tolist(),
            "radial": w.radial(k, channel.base.kappa).tolist(),
        })
    return waves


def tagged_model_dump(channel):
    model = tagged.TaggedModel(channel)  # defaults: k_max=1.2, nk=280, nc=96
    ms_ion = spin.m_values(channel.j_ion).tolist()
    k_pts = np.array([0.05, 0.10, 0.20, 0.30, 0.45, 0.60])
    c_pts = np.array([-0.9, -0.5, 0.0, 0.5, 0.9])

    # replicate TaggedModel's own nearest-below-cell lookup (tagged.py:256-258)
    ik = np.clip(np.searchsorted(model.k, k_pts) - 1, 0, model.k.size - 2)
    ic = np.clip(np.searchsorted(model.c, c_pts) - 1, 0, model.c.size - 2)

    n_entries = []
    struck_entries = []
    p2_entries = []
    for M in ms_ion:
        n_grid = model.n_of_kc(M)
        n_at_pts = [[float(n_grid[i, cc]) for cc in ic] for i in ik]
        n_entries.append({"M": M, "n_at_k_c": n_at_pts})

        struck = model.struck_populations(M)  # (n_mS, nk, nc)
        ms_struck = model._ms_struck.tolist()
        struck_at_pts = [[[float(struck[s_i, i, cc]) for cc in ic]
                          for i in ik] for s_i in range(len(ms_struck))]
        struck_entries.append({"M": M, "m_s_order": ms_struck,
                               "p_at_k_c": struck_at_pts})

        p2_entries.append({"M": M, "p2_moment": model.p2_moment(M)})

    out = {
        "grid": {"k_max": model.k[-1] * 0 + 1.2, "nk": model.k.size,
                "nc": model.c.size, "k_first": float(model.k[0]),
                "k_last": float(model.k[-1]), "c_first": float(model.c[0]),
                "c_last": float(model.c[-1])},
        "k_pts": k_pts.tolist(), "c_pts": c_pts.tolist(),
        "n_of_kc": n_entries,
        "struck_populations": struck_entries,
        "population_integrated": [
            {"M": M, "p_m_s": model.population_integrated(M).tolist()}
            for M in ms_ion],
        "norm": [{"M": M, "value": model.norm(M)} for M in ms_ion],
        "vector_dilution": model.vector_dilution(),
        "p2_moment": p2_entries,
        "p2_moment_mixture_uniform": model.p2_moment_mixture(
            [1.0 / len(ms_ion)] * len(ms_ion)),
    }
    if abs(channel.s_channel - 1.0) < 1e-9:
        out["tensor_dilution"] = model.tensor_dilution()
    else:
        MISSING.append(("TaggedModel.tensor_dilution channel=%s" %
                        channel.label,
                        "s_channel=%g != 1: tensor_dilution raises "
                        "ValueError by design (tagged.py:299-302)" %
                        channel.s_channel))
    return out


def boost_spectator_dump(channel, configs):
    kcp = [(0.05, 0.3, 0.0), (0.15, -0.5, 1.0), (0.30, 0.0, 3.0),
           (0.45, 0.7, 5.0)]
    ms_ion = spin.m_values(channel.j_ion).tolist()
    rows = []
    for cfg in configs:
        for k, c, phi_k in kcp:
            for theta_s, phi_s in ((0.0, 0.0), (0.5, 0.0)):
                lab = tagged.boost_spectator(
                    channel, k, c, phi_k, cfg["ion_momentum_per_nucleon"],
                    theta_s=theta_s, phi_s=phi_s)
                rows.append({
                    "electron_energy": cfg["electron_energy"],
                    "ion_momentum_per_nucleon":
                        cfg["ion_momentum_per_nucleon"],
                    "k": k, "c": c, "phi_k": phi_k,
                    "theta_s": theta_s, "phi_s": phi_s,
                    "pT": float(lab["pT"]), "theta": float(lab["theta"]),
                    "p_lab": float(lab["p_lab"]), "R": float(lab["R"]),
                    "xL": float(lab["xL"]), "kx": float(lab["kx"]),
                    "ky": float(lab["ky"]), "kz": float(lab["kz"]),
                    "phi_spec": float(lab["phi_spec"]),
                })
    return rows


PROVENANCE_NOTE_TAGGED = (
    "The li6_alpha and deuteron `model` blocks are dumped from the FIXED "
    "LiPolGen C++ library (validation/repin_tagged_from_lipolgen.py), not "
    "from polligen.  polligen's tagged._amp2_table omits the i^L "
    "partial-wave phase and therefore carries the INVERTED S-D interference "
    "sign; regenerating those blocks from it would re-bake the bug the "
    "2026-09-06 fix removed.  See "
    "docs/benchmarking/07_cw_sign_investigation.md and "
    "docs/open_items/run_2026-09-06/phase_CW_numbers.md.  Everything else in "
    "this file -- including li7_alpha's model block, one L = 1 wave with "
    "nothing to interfere with -- is still polligen's."
)


# Spin-1 tagged channels whose `model` block this script MUST NOT write.
#
# 2026-09-06: the tagged sector's S-D interference sign was inverted -- the
# partial-wave sum needs phi_L = i^L psi_L and had no phase at all.  The C++
# was fixed (src/core/tagged.cpp `build_amp2`); `polligen.tagged._amp2_table`
# (tagged.py:243-248) was NOT, and still carries the missing phase, so
# anything this script computes for an S+D channel's `model` block is the
# refuted sign.  Those two blocks are re-pinned from the FIXED LiPolGen by
# `validation/repin_tagged_from_lipolgen.py` and CARRIED THROUGH here.
#
# `li7_alpha` has one L = 1 wave -- the common i is a global phase, nothing
# interferes, nothing moved -- so its block is still polligen's and is
# rewritten normally.
#
# See docs/benchmarking/07_cw_sign_investigation.md.
TAGGED_MODEL_NOT_FROM_POLLIGEN = ("li6_alpha", "deuteron")


def _carried_tagged_models():
    """The re-pinned `model` blocks already on disk, keyed by channel."""
    path = OUT / "tagged.json"
    if not path.is_file():
        raise SystemExit(
            "validation/reference/tagged.json is missing and its %s model "
            "blocks cannot be regenerated from polligen (inverted S-D sign, "
            "see docs/benchmarking/07_cw_sign_investigation.md).  Restore the "
            "file, or rebuild it with "
            "validation/repin_tagged_from_lipolgen.py."
            % ", ".join(TAGGED_MODEL_NOT_FROM_POLLIGEN))
    with open(path) as fh:
        doc = json.load(fh)
    return {k: doc["channels"][k]["model"]
            for k in TAGGED_MODEL_NOT_FROM_POLLIGEN}


def build_tagged():
    carried = _carried_tagged_models()
    channels = {
        "li6_alpha": tagged.li6_alpha_channel(),
        "li7_alpha": tagged.li7_alpha_channel(),
        "deuteron": tagged.deuteron_channel(),
    }
    beam_species = {"li6_alpha": "6Li", "li7_alpha": "7Li", "deuteron": "d"}

    out = {"channels": {}, "P_D_LI6": tagged.P_D_LI6,
           "P_D_DEUTERON": tagged.P_D_DEUTERON,
           "provenance": "LiPolGen post-fix, formerly polligen",
           "provenance_note": PROVENANCE_NOTE_TAGGED}
    for key, channel in channels.items():
        cfgs = beams.default_configs(beam_species[key])
        cfg_dicts = [{"electron_energy": c.electron_energy,
                     "ion_momentum_per_nucleon": c.ion_momentum_per_nucleon}
                    for c in cfgs]
        out["channels"][key] = {
            "label": channel.label, "j_ion": channel.j_ion,
            "s_struck": channel.s_struck, "s_spec": channel.s_spec,
            "s_channel": channel.s_channel,
            "base": {
                "beam_A": channel.base.beam_A, "beam_Z": channel.base.beam_Z,
                "spectator": channel.base.spectator,
                "spectator_A": channel.base.spectator_A,
                "spectator_Z": channel.base.spectator_Z,
                "separation_energy": channel.base.separation_energy,
                "l_wave": channel.base.l_wave,
                "m_spec": channel.base.m_spec, "m_beam": channel.base.m_beam,
                "m_partner": channel.base.m_partner,
                "kappa": channel.base.kappa,
            },
            "waves": wave_dump(key, channel),
            "beam_configs": cfg_dicts,
            "model": (carried[key] if key in carried
                      else tagged_model_dump(channel)),
            "boost_spectator": boost_spectator_dump(channel, cfg_dicts),
        }
    return out


doc("tagged.json",
    "`polligen.tagged` (PolarizedLithiumSim/evgen/polligen/tagged.py) using\n"
    "the DEFAULT channel constructors (all default beta=0.30, and for 6Li\n"
    "p_d=tagged.P_D_LI6, for the deuteron control p_d=tagged.P_D_DEUTERON).\n\n"
    "**THIS FILE NO LONGER TRACKS polligen FOR THE SPIN-1 MODEL BLOCKS "
    "(2026-09-06).**\n"
    "`polligen.tagged._amp2_table` (tagged.py:243-248) sums the partial "
    "waves with\n"
    "NO `i^L`: it feeds `psi_L` into an amplitude that needs "
    "`phi_L = i^L psi_L`.\n"
    "For an S+D channel that is not a global phase -- it is +1 on L=0 and -1 "
    "on\n"
    "L=2 -- so polligen's tagged sector carries the S-D interference sign\n"
    "INVERTED, against Cosyn-Weiss II Eq. (6.12) (by up to 2.74 in an "
    "asymmetry\n"
    "whose whole range is [-2, 1]), against LiPolGen's own deuteron "
    "quadrupole\n"
    "sign gate, and against LiPolGen's own b1 sector, which applies the "
    "phase\n"
    "explicitly.  LiPolGen fixed it on 2026-09-06 "
    "(`src/core/tagged.cpp`\n"
    "`build_amp2`, one `(-1)^floor(L/2)`); polligen was NOT touched and "
    "still\n"
    "carries the bug.  Regenerating `channels.li6_alpha.model` or\n"
    "`channels.deuteron.model` from polligen would therefore re-bake the "
    "refuted\n"
    "sign and the rtol-1e-12 gate would go on certifying it.\n\n"
    "Those two blocks are instead RE-PINNED from the fixed C++ library by\n"
    "`validation/repin_tagged_from_lipolgen.py` (provenance `\"LiPolGen post-fix,\n"
    "formerly polligen\"`, recorded in the file's own `provenance` /\n"
    "`provenance_note` keys), and `dump_polligen_reference.py` carries them\n"
    "through unchanged rather than overwriting them.  What moved in the re-pin,\n"
    "MEASURED 2026-09-06 as the re-pinned file against the pre-fix dump\n"
    "(`git show HEAD:validation/reference/tagged.json`), not copied from the\n"
    "investigation's own table: `n_of_kc` (up to +725% on 6Li, +19215% on the\n"
    "AV18 deuteron control), `struck_populations` (up to 0.86 absolute),\n"
    "`p2_moment` (sign flip, x1.28 to x2.03), `p2_moment_mixture_uniform`\n"
    "(-5.3160743e-05 -> -5.2153460e-05 on 6Li, 1.9e-2 rel; -5.3694953e-05 ->\n"
    "-5.3337762e-05 on the deuteron, 6.7e-3 rel -- it is gated at 1e-9, so it had\n"
    "to be re-pinned too), `norm` (up to 5.8e-5 rel: 6Li M=0 1.000026689626 ->\n"
    "0.999968571520; 4.0e-5 on the deuteron), `population_integrated` (up to\n"
    "3.0e-6 abs: 6Li M=0 0.947966373524 -> 0.947963349333) and the dilutions\n"
    "(<=4.3e-6 rel).  Until 2026-09-06 this list read `norm` \"+2e-5\" and\n"
    "`population_integrated` \"<=5e-7 abs\" and did not mention\n"
    "`p2_moment_mixture_uniform` at all -- three figures taken from\n"
    "`07_cw_sign_investigation.md` section 6.1 rather than measured here, where\n"
    "that section records the last one as not moving.  Reason, derivations and\n"
    "the full before/after tables: `docs/benchmarking/07_cw_sign_investigation.md`\n"
    "and `docs/open_items/run_2026-09-06/phase_CW_numbers.md`.\n\n"
    "EVERYTHING ELSE IN THE FILE IS STILL polligen's, untouched: `waves`,\n"
    "`base`, `beam_configs`, `boost_spectator`, `P_D_LI6`, `P_D_DEUTERON`, "
    "the\n"
    "channel scalars, and the WHOLE of `channels.li7_alpha.model` -- 7Li "
    "alpha-tag\n"
    "is a single L=1 wave, so the common `i` is a global phase, and the fix\n"
    "moves its `n_of_kc` and `p2_moment` by exactly zero (measured, "
    "bit for bit).\n"
    "The re-pin script re-checks that block against the live C++ at rtol "
    "1e-12\n"
    "instead of overwriting it.  No other reference JSON moved.\n\n"
    "`channels.<li6_alpha|li7_alpha|deuteron>`: built by\n"
    "`tagged.li6_alpha_channel()` / `li7_alpha_channel()` / "
    "`deuteron_channel()` (tagged.py:181,188,195). `base` is the underlying\n"
    "`polli_fastsim.spectator.ClusterChannel` (spectator.py:111): `m_spec`\n"
    "(spectator.py:124), `m_beam` (spectator.py:128, physical nuclear mass),\n"
    "`m_partner` (spectator.py:143), `kappa` (spectator.py:168, "
    "sqrt(2*mu*S)).\n\n"
    "`waves[i]`: one `tagged.Wave` (tagged.py:126) of `channel.waves` --\n"
    "`l`, `prob`, `beta`, and `radial` = `Wave.radial(k_grid, kappa)`\n"
    "(tagged.py:140; L=0 Hulthen 1/(k^2+kappa^2)-1/(k^2+beta^2), L=1\n"
    "k/((k^2+kappa^2)(k^2+beta^2)), L=2 k^2/((k^2+kappa^2)(k^2+beta^2)^2)) on\n"
    "`k_grid` = 40 points linspace(1e-4, 1.2, 40) GeV -- UNNORMALIZED (the\n"
    "TaggedModel normalizes internally, see below).\n\n"
    "`model`: `tagged.TaggedModel(channel)` (tagged.py:203) at its DEFAULT\n"
    "grid (k_max=1.2, nk=280, nc=96 -- `grid` records the exact edges/sizes,\n"
    "which the C++ port must reproduce bit-for-bit since `n_of_kc` looks up\n"
    "the NEAREST cell at or below the query point, `np.clip(np.searchsorted"
    "(model.k, k)-1, 0, nk-2)`, tagged.py:256-258 -- this script replicates\n"
    "that exact lookup in Python to pick `k_pts`/`c_pts` cell values, so the\n"
    "reference values are themselves grid CELL values, not interpolated).\n"
    "  * `n_of_kc[i].n_at_k_c`: `TaggedModel.n_of_kc(M)` (tagged.py:248) at\n"
    "    `k_pts` x `c_pts` (6x5), indexed [k][c].\n"
    "  * `struck_populations[i].p_at_k_c`: `TaggedModel.struck_populations(M)`\n"
    "    (tagged.py:260), shape (n_mS, nk, nc), sliced at the same\n"
    "    (k_pts, c_pts) cells; `m_s_order` gives the m_S ordering (+S_c...-S_c).\n"
    "  * `population_integrated[i].p_m_s`: `TaggedModel.population_integrated(M)`\n"
    "    (tagged.py:266), the k/khat-integrated channel-spin populations.\n"
    "  * `norm[i].value`: `TaggedModel.norm(M)` (tagged.py:272), grid\n"
    "    quadrature of the normalization integral (should be ~1, NOT exact\n"
    "    -- it is itself a finite-grid quadrature, so match this one at a\n"
    "    looser tolerance, e.g. 1e-6, not rtol=1e-12).\n"
    "  * `vector_dilution`: `TaggedModel.vector_dilution()` (tagged.py:292).\n"
    "  * `tensor_dilution` (only present when s_channel=1, i.e. li6_alpha and\n"
    "    deuteron -- s_channel=0.5 for li7_alpha, where `tensor_dilution`\n"
    "    raises `ValueError` by design, tagged.py:299-302; see the JSON's\n"
    "    top-level `_could_not_call` list): `TaggedModel.tensor_dilution()`\n"
    "    (tagged.py:299).\n"
    "  * `p2_moment[i]`: `TaggedModel.p2_moment(M)` (tagged.py:307), per M.\n"
    "  * `p2_moment_mixture_uniform`: `TaggedModel.p2_moment_mixture(pops)`\n"
    "    (tagged.py:314) at a uniform fill (also a grid quadrature; same\n"
    "    looser-tolerance caveat as `norm`).\n\n"
    "`boost_spectator`: `tagged.boost_spectator(channel, k, c, phi_k,\n"
    "p_per_nucleon, theta_s, phi_s)` (tagged.py:335), an EXACT closed-form\n"
    "boost (no grid quadrature -- rtol=1e-12 applies), at a fixed set of\n"
    "(k, c, phi_k) x (theta_s, phi_s) x each of the channel's 3\n"
    "`beam_configs` (`beams.default_configs` for the channel's own beam\n"
    "species: 6Li, 7Li, or d). Fields: pT, theta, p_lab, R (rigidity ratio),\n"
    "xL, kx, ky, kz (spin-frame components after rotation), phi_spec (lab\n"
    "azimuth).")


# ===========================================================================
# 5. spectator.json
# ===========================================================================

def build_spectator():
    named_channels = {
        "DEUTERON_P_TAG": spectator.DEUTERON_P_TAG,
        "DEUTERON_N_TAG": spectator.DEUTERON_N_TAG,
        "HE3_P_TAG": spectator.HE3_P_TAG,
        "LI6_ALPHA_TAG": spectator.LI6_ALPHA_TAG,
        "LI6_D_TAG": spectator.LI6_D_TAG,
        "LI7_ALPHA_TAG": spectator.LI7_ALPHA_TAG,
        "LI7_T_TAG": spectator.LI7_T_TAG,
    }
    k_grid = np.linspace(1e-4, 1.2, 40)
    p_per_nucleon_list = [40.8, 99.5, 137.5]  # representative 6Li/7Li values
    kxyz_list = [(0.0, 0.0, 0.0), (0.1, 0.0, 0.0), (0.0, 0.1, 0.2),
                (-0.05, 0.07, -0.03)]

    out = {"constants": {
        "M_U": spectator.M_U,
        "MASSES": dict(spectator.MASSES),
        "NUCLEUS_MASS": {"%d_%d" % k: v
                         for k, v in spectator.NUCLEUS_MASS.items()},
    }, "channels": {}}

    for name, ch in named_channels.items():
        r_at_zero_formula = (float("nan") if ch.spectator_Z == 0 else
                             (ch.m_spec / ch.spectator_Z)
                             / (ch.m_beam / ch.beam_Z))
        boosts = []
        for p_u in p_per_nucleon_list:
            for kx, ky, kz in kxyz_list:
                res = spectator._boost_fragment(
                    ch, p_u, np.array([kx]), np.array([ky]), np.array([kz]),
                    ch.m_spec, ch.spectator_Z, ch.spectator_A)
                boosts.append({
                    "p_per_nucleon": p_u, "kx": kx, "ky": ky, "kz": kz,
                    "pT": float(res["pT"][0]), "theta": float(res["theta"][0]),
                    "phi": float(res["phi"][0]), "p_lab": float(res["p_lab"][0]),
                    "R": float(res["R"][0]), "xL": float(res["xL"][0]),
                    "k": float(res["k"][0]),
                })
        out["channels"][name] = {
            "beam_A": ch.beam_A, "beam_Z": ch.beam_Z,
            "spectator": ch.spectator, "spectator_A": ch.spectator_A,
            "spectator_Z": ch.spectator_Z,
            "separation_energy": ch.separation_energy, "l_wave": ch.l_wave,
            "m_spec": ch.m_spec, "m_beam": ch.m_beam,
            "m_partner": ch.m_partner, "kappa": ch.kappa,
            "R_at_k_zero": r_at_zero_formula,
            "momentum_density": {
                "beta": 0.30, "k_grid": k_grid.tolist(),
                "n_of_k": spectator.momentum_density(
                    k_grid, ch.kappa, 0.30, ch.l_wave).tolist(),
            },
            "boost_fragment": boosts,
        }
    return out


doc("spectator.json",
    "`polli_fastsim.spectator` "
    "(PolarizedLithiumSim/fastsim/polli_fastsim/spectator.py).\n\n"
    "`constants`: `M_U` (spectator.py:43), `MASSES` (spectator.py:44, keyed\n"
    "by particle name), `NUCLEUS_MASS` (spectator.py:75, keyed `\"<Z>_<A>\"`,\n"
    "the AME2020-derived physical nuclear masses).\n\n"
    "`channels.<name>`: one of `spectator.{DEUTERON_P_TAG, DEUTERON_N_TAG,\n"
    "HE3_P_TAG, LI6_ALPHA_TAG, LI6_D_TAG, LI7_ALPHA_TAG, LI7_T_TAG}`\n"
    "(spectator.py:187-207), each a `ClusterChannel` (spectator.py:111):\n"
    "m_spec (spectator.py:124), m_beam (spectator.py:128), m_partner\n"
    "(spectator.py:143), kappa (spectator.py:168, sqrt(2*mu*S) with mu the\n"
    "reduced mass of the FREE spectator+partner). `R_at_k_zero` is the\n"
    "closed form `(m_spec/spectator_Z)/(m_beam/beam_Z)` (spectator.py\n"
    "module docstring, verified against `_boost_fragment` at k=0 --\n"
    "independent of p_per_nucleon, see `boost_fragment` rows with\n"
    "kx=ky=kz=0).\n\n"
    "`momentum_density.n_of_k`: `spectator.momentum_density(k, kappa,\n"
    "beta=0.30, l_wave)` (spectator.py:212, the UNNORMALIZED |psi(k)|^2 --\n"
    "Hulthen^2 for l_wave=0, [k/((k^2+kappa^2)(k^2+beta^2))]^2 for l_wave=1)\n"
    "on `k_grid` (40 points, linspace(1e-4, 1.2, 40) GeV).\n\n"
    "`boost_fragment`: `spectator._boost_fragment(channel, p_per_nucleon,\n"
    "kx, ky, kz, m=channel.m_spec, frag_Z=channel.spectator_Z,\n"
    "frag_A=channel.spectator_A)` (spectator.py:238) -- i.e. exactly what\n"
    "`spectator_lab_kinematics` (spectator.py:280) calls internally, but at\n"
    "FIXED (kx, ky, kz) rather than a random `sample_k` draw: reproducing\n"
    "`spectator_lab_kinematics`'s own random sampling bit-for-bit would\n"
    "require replicating numpy's PCG64 bit generator in C++, which is out\n"
    "of scope for a physics cross-check, so this table isolates the\n"
    "boost/kinematics formula itself (the only thing `sample_k` feeds it).\n"
    "Evaluated at p_per_nucleon in {40.8, 99.5, 137.5} GeV (representative\n"
    "6Li/7Li values) x 4 fixed (kx, ky, kz) [GeV] triples including (0,0,0).\n"
    "Fields: pT, theta, phi (lab azimuth), p_lab, R (rigidity ratio), xL, k.")


# ===========================================================================
# 6. coherent.json
# ===========================================================================

def build_coherent():
    default = coherent.CoherentScenario()
    out = {
        "constants": {
            "M_LI6": coherent.M_LI6, "GEV_PER_FM_INV": coherent.GEV_PER_FM_INV,
            "RATE_WEIGHT_SYST": coherent.RATE_WEIGHT_SYST,
            "MANTYSAARI_A2_DEUTERON": {str(k): list(v) for k, v in
                                      coherent.MANTYSAARI_A2_DEUTERON.items()},
        },
        "scenario_defaults": {
            "f0": default.f0, "x_coh": default.x_coh,
            "slope_b": default.slope_b, "amp": default.amp,
            "eps_b0": default.eps_b0,
        },
        "gaussian_slope": [
            {"r_rms_fm": r, "B": coherent.gaussian_slope(r)}
            for r in (2.32, 2.45, 2.589)
        ],
    }

    xs = [1e-4, 1e-3, 0.01, 0.02, 0.05, 0.1, 0.3]
    out["coherent_fraction"] = [
        {"x": x, "value": float(default.coherent_fraction(x))} for x in xs
    ]

    pt_cuts = [0.1, 0.2, 0.3, 0.45]
    out["tag_acceptance"] = [
        {"pt_cut": p, "value": default.tag_acceptance(p),
         "mean_t_tagged": default.mean_t_tagged(p)} for p in pt_cuts
    ]

    sigma_thetas = [73e-6, 164e-6, 220e-6, 380e-6]
    p_per_nucleons = [40.8, 99.5, 137.5]
    out["tag_acceptance_angular"] = [
        {"sigma_theta": st, "p_per_nucleon": pu, "a_beam": 6,
         "n_sigma": 10.0,
         "value": float(default.tag_acceptance_angular(st, pu))}
        for st in sigma_thetas for pu in p_per_nucleons
    ]

    t_abs_vals = [0.0, 0.02, 0.05, 0.1, 0.2]
    pzz_vals = [0.6, -0.6, 1.0]
    out["a2_deformation"] = [
        {"t_abs": t, "pzz": pzz,
         "a2_deformation": float(default.a2_deformation(t, pzz)),
         "cos2phi_coefficient_deformation":
             float(default.cos2phi_coefficient_deformation(t, pzz))}
        for t in t_abs_vals for pzz in pzz_vals
    ]
    out["a2_tagged"] = [
        {"pt_cut": p, "pzz": pzz, "value": default.a2_tagged(p, pzz)}
        for p in pt_cuts for pzz in pzz_vals
    ]

    t_list = [0.0, 0.01, 0.05, 0.1, 0.3]
    phi_t_list = [0.0, 1.0, 3.0]
    p_u = 99.5
    recoil_rows = []
    for t_abs in t_list:
        for phi_t in phi_t_list:
            r = coherent.recoil_lab(t_abs, phi_t, p_u)
            recoil_rows.append({
                "t_abs": t_abs, "phi_t": phi_t, "p_per_nucleon": p_u,
                "pT": float(r["pT"]), "theta": float(r["theta"]),
                "R": float(r["R"]), "xL": float(r["xL"]),
            })
    out["recoil_lab"] = recoil_rows

    frag_rows = []
    for beam_a, beam_z, table in ((6, 3, coherent.LI6_BREAKUP),
                                  (7, 3, coherent.LI7_BREAKUP)):
        for name, thr, frags in table:
            for fname, a, z in frags:
                frag_rows.append({
                    "beam_a": beam_a, "beam_z": beam_z, "channel": name,
                    "threshold_mev": thr, "fragment": fname, "a": a, "z": z,
                    "rigidity": coherent.fragment_rigidity(
                        a, z, beam_a=beam_a, beam_z=beam_z),
                })
    out["fragment_rigidity"] = frag_rows
    out["breakup_tables"] = {
        "6Li": [{"name": n, "threshold_mev": t,
                "fragments": [{"name": fn, "a": a, "z": z}
                             for fn, a, z in fr]}
               for n, t, fr in coherent.LI6_BREAKUP],
        "7Li": [{"name": n, "threshold_mev": t,
                "fragments": [{"name": fn, "a": a, "z": z}
                             for fn, a, z in fr]}
               for n, t, fr in coherent.LI7_BREAKUP],
    }
    return out


doc("coherent.json",
    "`polligen.coherent` (PolarizedLithiumSim/evgen/polligen/coherent.py).\n\n"
    "`constants`: `M_LI6` (coherent.py:47), `GEV_PER_FM_INV` (coherent.py:48,\n"
    "hbar*c), `RATE_WEIGHT_SYST` (coherent.py:213), `MANTYSAARI_A2_DEUTERON`\n"
    "(coherent.py:221, keyed by |t| [GeV^2] as a string, value = (a2 m=0,\n"
    "a2 m=+-1)).\n\n"
    "`scenario_defaults`: the field defaults of `coherent.CoherentScenario`\n"
    "(coherent.py:57): f0=0.04, x_coh=0.01, slope_b=50.0, amp=0.01,\n"
    "eps_b0=-0.08.\n\n"
    "`gaussian_slope`: `coherent.gaussian_slope(r_rms_fm)` (coherent.py:51,\n"
    "= (r_rms_fm/GEV_PER_FM_INV)^2 / 3) at the matter/charge radii quoted in\n"
    "the docstring.\n\n"
    "`coherent_fraction`: `CoherentScenario.coherent_fraction(x)`\n"
    "(coherent.py:111, = f0/(1+(x/x_coh)^2)) at the default scenario.\n\n"
    "`tag_acceptance`: `CoherentScenario.tag_acceptance(pt_cut)`\n"
    "(coherent.py:123, = exp(-B*pt_cut^2)) and `mean_t_tagged(pt_cut)`\n"
    "(coherent.py:128, = pt_cut^2 + 1/B).\n\n"
    "`tag_acceptance_angular`: `CoherentScenario.tag_acceptance_angular(\n"
    "sigma_theta, p_per_nucleon, a_beam=6, n_sigma=10.0)` (coherent.py:132,\n"
    "= exp(-B*(n_sigma*sigma_theta*a_beam*p_per_nucleon)^2)).\n\n"
    "`a2_deformation`: `CoherentScenario.a2_deformation(t_abs, pzz)`\n"
    "(coherent.py:155, = -(pzz/4)*eps_b0*B*|t|) and\n"
    "`cos2phi_coefficient_deformation(t_abs, pzz)` (coherent.py:178,\n"
    "= 2*a2_deformation).\n\n"
    "`a2_tagged`: `CoherentScenario.a2_tagged(pt_cut, pzz)` (coherent.py:196,\n"
    "= a2_deformation(mean_t_tagged(pt_cut), pzz)).\n\n"
    "`recoil_lab`: `coherent.recoil_lab(t_abs, phi_t, p_per_nucleon,\n"
    "x_pom=0.0)` (coherent.py:229, exact closed-form kinematics of the\n"
    "intact-6Li recoil): pT, theta, R, xL.\n\n"
    "`fragment_rigidity`: `coherent.fragment_rigidity(a, z, beam_a, beam_z)`\n"
    "(coherent.py:352, = (nucleus_mass(z,a)/z) / (nucleus_mass(beam_z,\n"
    "beam_a)/beam_z)) for every fragment of `coherent.LI6_BREAKUP` /\n"
    "`LI7_BREAKUP` (coherent.py:306,331), whose raw contents are also given\n"
    "verbatim in `breakup_tables`.")


# ===========================================================================
# 7. bookkeeping.json
# ===========================================================================

def category_dump(cat):
    out = {
        "name": cat.name, "j": cat.j, "populations": list(cat.populations),
        "lam_e": cat.lam_e, "pe": cat.pe, "theta_s": cat.theta_s,
        "phi_s": cat.phi_s, "lumi_fraction": cat.lumi_fraction,
    }
    try:
        out["moments"] = list(cat.moments())
    except ValueError:
        # spin.moments_along_axis (spin.py:165) is defined for j=1, 3/2
        # only; for j=1/2 the vector moment alone is well defined (no
        # rank-2 for a two-level system), computed directly here.
        ms = spin.m_values(cat.j)
        pops = np.asarray(cat.populations, dtype=float)
        out["moments"] = [float((ms * pops).sum() / cat.j)]
        MISSING.append(("SpinCategory.moments() name=%s (j=%g)" %
                        (cat.name, cat.j),
                        "spin.moments_along_axis raises ValueError for "
                        "j not in {1, 1.5} (spin.py:177); vector moment "
                        "computed directly instead, tensor/octupole n/a"))
    return out


def plan_dump(plan, total_lumi_pb=1e6):
    return {
        "categories": [category_dump(c) for c in plan.categories],
        "pe_true": plan.pe_true, "pz_true": plan.pz_true,
        "pzz_true": plan.pzz_true,
        "lumi_shares_at_1e6_pb": plan.lumi_shares(total_lumi_pb),
    }


def build_bookkeeping():
    out = {"plans": {}, "rel_lumi_bias": {}}

    out["plans"]["helicity_flip_j1_maxent_anchor"] = plan_dump(
        bookkeeping.helicity_flip_plan(1.0, 8.0 / 13.0, PE))
    out["plans"]["helicity_flip_j32_maxent_anchor"] = plan_dump(
        bookkeeping.helicity_flip_plan(1.5, 0.7, PE))
    out["plans"]["helicity_flip_j12"] = plan_dump(
        bookkeeping.helicity_flip_plan(0.5, 0.7, PE))
    out["plans"]["helicity_flip_j1_explicit_pzz"] = plan_dump(
        bookkeeping.helicity_flip_plan(1.0, 0.6, PE, pzz=0.35))
    out["plans"]["helicity_flip_j1_tilted_offset"] = plan_dump(
        bookkeeping.helicity_flip_plan(1.0, 0.7, PE, theta_s=0.9,
                                       phi_s=0.4, rel_lumi_offset=1e-4))
    out["plans"]["tensor_thirds_j1"] = plan_dump(
        bookkeeping.tensor_thirds_plan(0.7, 0.6, rel_lumi_offset=1e-4))
    out["plans"]["transverse_tensor_j1"] = plan_dump(
        bookkeeping.transverse_tensor_plan(0.6))
    out["plans"]["tensor_flip_j1"] = plan_dump(
        bookkeeping.tensor_flip_plan(0.6, share_plus=0.5,
                                     rel_lumi_offset=1e-4))

    with_off = bookkeeping.with_offset(
        bookkeeping.tensor_thirds_plan(0.7, 0.6), "azz0", 1e-4)
    out["plans"]["tensor_thirds_j1_with_offset_azz0"] = plan_dump(with_off)

    # spin-temperature ladder rational anchors (t = 3), from test_bookkeeping.py
    out["spin_temperature_ladder_t3"] = {
        "j1": {"t": 3.0, "populations": [9.0 / 13.0, 3.0 / 13.0, 1.0 / 13.0],
              "pz": 8.0 / 13.0, "pzz": 4.0 / 13.0},
        "j32": {"t": 3.0,
               "populations": [27.0 / 40.0, 9.0 / 40.0, 3.0 / 40.0, 1.0 / 40.0],
               "pz": 0.7, "t_moment": 0.4},
    }

    deltas = [1e-4]
    pzz_vals = [0.6, 0.4]
    pe_vals = [0.7]
    pz_vals = [0.7, 0.6]
    out["rel_lumi_bias"]["azz_rel_lumi_bias"] = [
        {"offset": d, "pzz": pzz,
         "value": bookkeeping.azz_rel_lumi_bias(d, pzz)}
        for d in deltas for pzz in pzz_vals
    ]
    out["rel_lumi_bias"]["apar_rel_lumi_bias"] = [
        {"offset": d, "pe": pe, "pz": pz,
         "value": bookkeeping.apar_rel_lumi_bias(d, pe, pz)}
        for d in deltas for pe in pe_vals for pz in pz_vals
    ]
    return out


doc("bookkeeping.json",
    "`polligen.bookkeeping` "
    "(PolarizedLithiumSim/evgen/polligen/bookkeeping.py).\n\n"
    "`plans.<name>`: the four standard run-plan constructors, each a list of\n"
    "`SpinCategory` (bookkeeping.py:36: name, j, populations, lam_e, pe,\n"
    "theta_s, phi_s, lumi_fraction, and `moments()` = "
    "`spin.moments_along_axis(j, populations)`) plus the `RunPlan`\n"
    "(bookkeeping.py:52) true polarizations and `lumi_shares(1e6)`\n"
    "(bookkeeping.py:74, absolute lumi at an arbitrary reference total of\n"
    "1e6 pb^-1 = 1 fb^-1 -- purely `lumi_fraction * total`, exact rational\n"
    "arithmetic).\n"
    "  * `helicity_flip_j1_maxent_anchor` /\n"
    "    `helicity_flip_j32_maxent_anchor`: `bookkeeping.helicity_flip_plan(\n"
    "    j, pz, pe)` (bookkeeping.py:80) with pzz=None (spin-temperature\n"
    "    fill via `spin.populations_maxent`) at the t=3 rational anchors\n"
    "    pz=8/13 (J=1, expect pzz_true=4/13) and pz=0.7 (J=3/2, expect\n"
    "    pzz_true=0.4) -- see `spin_temperature_ladder_t3` below for the\n"
    "    exact closed-form populations (bisection-derived in the module, so\n"
    "    match this ladder at a looser tolerance than 1e-12, e.g. 1e-10).\n"
    "  * `helicity_flip_j12`: J=1/2 branch (no rank-2 moment, populations\n"
    "    (1+pz)/2, (1-pz)/2 exactly).\n"
    "  * `helicity_flip_j1_explicit_pzz`: pzz=0.35 explicit (closed-form\n"
    "    `spin.spin1_populations`, exact to double precision).\n"
    "  * `helicity_flip_j1_tilted_offset`: theta_s=0.9, phi_s=0.4,\n"
    "    rel_lumi_offset=1e-4 (boosts the lam_e=+1 share only, see\n"
    "    bookkeeping.py:102).\n"
    "  * `tensor_thirds_j1`: `bookkeeping.tensor_thirds_plan(pz=0.7,\n"
    "    pzz=0.6, rel_lumi_offset=1e-4)` (bookkeeping.py:114): 3 categories\n"
    "    (+pzz, -pzz, m0-enriched -2pzz).\n"
    "  * `transverse_tensor_j1`: `bookkeeping.transverse_tensor_plan(\n"
    "    pzz=0.6)` (bookkeeping.py:140): 1 category at theta_S=pi/2.\n"
    "  * `tensor_flip_j1`: `bookkeeping.tensor_flip_plan(pzz=0.6,\n"
    "    share_plus=0.5, rel_lumi_offset=1e-4)` (bookkeeping.py:150).\n"
    "  * `tensor_thirds_j1_with_offset_azz0`: `bookkeeping.with_offset(\n"
    "    tensor_thirds_plan(0.7, 0.6), \"azz0\", 1e-4)` (bookkeeping.py:178) --\n"
    "    only the \"azz0\" category's lumi_fraction is scaled by (1+1e-4).\n\n"
    "`spin_temperature_ladder_t3`: the closed-form geometric-ladder\n"
    "populations p_m ~ t^m at t=3 for J=1 and J=3/2 (exact rationals, from\n"
    "the docstring/tests of bookkeeping.py and test_bookkeeping.py), the\n"
    "independent construction the maxent anchors above should match.\n\n"
    "`rel_lumi_bias`: `bookkeeping.azz_rel_lumi_bias(offset, pzz)`\n"
    "(bookkeeping.py:193, = -(2/3)*offset/pzz) and\n"
    "`bookkeeping.apar_rel_lumi_bias(offset, pe, pz)` (bookkeeping.py:199,\n"
    "= offset/(2*pe*pz)) at offset=1e-4 (the plans/05 reference delta).")


# ===========================================================================
# main
# ===========================================================================

def main():
    written = []

    written.append(write_json("spin.json", build_spin()))
    written.append(write_json("xsec.json", build_xsec()))
    written.append(write_json("beams.json", build_beams()))
    written.append(write_json("tagged.json", build_tagged()))
    written.append(write_json("spectator.json", build_spectator()))
    written.append(write_json("coherent.json", build_coherent()))
    written.append(write_json("bookkeeping.json", build_bookkeeping()))

    # write the "could not call" list into every file that had one, plus
    # print a summary; also emit it as its own small manifest.
    #
    # MERGE, do not overwrite: `validation/reference/` also holds files this
    # script does not make -- `b1_default_li6.json` comes from
    # `dump_b1_default_li6.py` -- and rewriting the manifest wholesale from
    # `written` silently dropped them (and the `generators` map that says
    # which script owns which file).  Only this script's own entries are
    # refreshed; foreign ones are carried through.
    manifest_path = OUT / "_manifest.json"
    doc_ = {}
    if manifest_path.is_file():
        with open(manifest_path) as fh:
            doc_ = json.load(fh)
    mine = [p.name for p in written]
    files = mine + [f for f in doc_.get("files", []) if f not in mine]
    gens = dict(doc_.get("generators", {}))
    for f in mine:
        # tagged.json's spin-1 model blocks are re-pinned from the fixed C++
        # (TAGGED_MODEL_NOT_FROM_POLLIGEN above) and only carried through
        # here, so the file's OWNER stays the re-pin script.  Claiming it
        # would send the next reader back to polligen's inverted S-D sign.
        if f == "tagged.json":
            gens.setdefault(f, "LiPolGen/validation/"
                               "repin_tagged_from_lipolgen.py")
            continue
        gens[f] = "LiPolGen/validation/dump_polligen_reference.py"
    write_json("_manifest.json", {
        "generator": "LiPolGen/validation/dump_polligen_reference.py",
        "polligen_version": polligen.__version__,
        "could_not_call": [{"what": w, "reason": r} for w, r in MISSING],
        "files": files,
        "generators": gens,
    })

    readme_path = HERE / "README.md"
    with open(readme_path, "w") as fh:
        fh.write("# polligen reference tables\n\n")
        fh.write(
            "Numeric reference tables dumped from the Python event "
            "generator `polligen`\n"
            "(`PolarizedLithiumSim/evgen/polligen`) by "
            "`dump_polligen_reference.py`, for a\n"
            "C++ reimplementation (LiPolGen) to be checked against at "
            "`rtol=1e-12`\n"
            "unless a table's description says otherwise (a few "
            "quantities in polligen\n"
            "are themselves iterative/quadrature results, e.g. "
            "`spin.populations_maxent`,\n"
            "`TaggedModel.norm`, `p2_moment_mixture` -- those are called "
            "out below).\n\n"
            "Every JSON file has an `\"inputs\"`-shaped section (grids, "
            "axes, kernel /\n"
            "scenario configuration, constants) and an outputs section "
            "with the actual\n"
            "computed numbers, organized per table below. Floats are "
            "written with Python's\n"
            "`repr()`-precision JSON serialization, i.e. full IEEE-754 "
            "double precision\n"
            "round-trips exactly.\n\n"
            "Regenerate with:\n\n"
            "```\n"
            "python3 LiPolGen/validation/dump_polligen_reference.py\n"
            "```\n\n"
            "This script only *imports* `polligen`/`polli_fastsim` from "
            "`PolarizedLithiumSim`;\n"
            "it does not modify anything there.\n\n"
        )
        for title, text in README_SECTIONS:
            fh.write("## %s\n\n%s\n\n" % (title, text))
        if MISSING:
            fh.write("## Functions this script could not call\n\n")
            for w, r in MISSING:
                fh.write("- `%s`: %s\n" % (w, r))
            fh.write("\n")

    sizes = [(p.name, p.stat().st_size) for p in sorted(written)]
    print("Wrote %d JSON files to %s:" % (len(written), OUT))
    for name, size in sizes:
        print("  %-20s %8d bytes" % (name, size))
    print("Wrote README to %s" % readme_path)
    if MISSING:
        print("\nCould not call:")
        for w, r in MISSING:
            print("  - %s: %s" % (w, r))
    else:
        print("\nEvery targeted function was called successfully.")


if __name__ == "__main__":
    main()
